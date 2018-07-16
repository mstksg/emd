{-# LANGUAGE BangPatterns                             #-}
{-# LANGUAGE GADTs                                    #-}
{-# LANGUAGE LambdaCase                               #-}
{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE TypeInType                               #-}
{-# LANGUAGE TypeOperators                            #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise       #-}

module Numeric.EMD (
    emd, emdTrace, emd'
  , EMD(..)
  , SiftCondition(..), defaultSC
  , sift
  -- * Debug
  , extrema
  , envelopes
  ) where

import           Control.Monad
import           Control.Monad.IO.Class
import           Data.Bifunctor
import           Data.Finite
import           Data.Functor.Identity
import           GHC.TypeNats
import           Numeric.Spline
import           Text.Printf
import qualified Data.Map                  as M
import qualified Data.Vector.Generic       as VG
import qualified Data.Vector.Generic.Sized as SVG

data SiftCondition a = SCStdDev a
                     | SCTimes Int
                     | SCOr (SiftCondition a) (SiftCondition a)
                     | SCAnd (SiftCondition a) (SiftCondition a)
  deriving Show

defaultSC :: Fractional a => SiftCondition a
defaultSC = SCStdDev 0.3

-- | 'True' if stop
testCondition
    :: (VG.Vector v a, Fractional a, Ord a)
    => SiftCondition a
    -> Int
    -> SVG.Vector v n a
    -> SVG.Vector v n a
    -> Bool
testCondition = \case
    SCStdDev t -> \_ v v' ->
      let sd = SVG.sum $ SVG.zipWith (\x x' -> (x-x')^(2::Int) / (x + eps)^(2::Int)) v v'
      in  sd <= t
    SCTimes l  -> \i _ _ -> i >= l
    SCOr  f g -> \i v v' -> testCondition f i v v' || testCondition g i v v'
    SCAnd f g -> \i v v' -> testCondition f i v v' && testCondition g i v v'
  where
    eps = 0.0000001

data EMD v n a = EMD { emdIMFs     :: ![SVG.Vector v n a]
                     , emdResidual :: !(SVG.Vector v n a)
                     }
  deriving Show

-- | EMD decomposition (Hilbert-Huang Transform) of a given time series
-- with a given sifting stop condition.
emd :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SiftCondition a
    -> SVG.Vector v (n + 1) a
    -> EMD v (n + 1) a
emd sc = runIdentity . emd' (const (pure ())) sc

-- | 'emd', but tracing results as IMFs are found.
emdTrace
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a, MonadIO m)
    => SiftCondition a
    -> SVG.Vector v (n + 1) a
    -> m (EMD v (n + 1) a)
emdTrace = emd' $ \case
    SRResidual _ -> liftIO $ putStrLn "Residual found!"
    SRIMF _ i    -> liftIO $ printf "IMF found (%d iterations)\n" i

-- | 'emd' with an optional callback for each found IMF.
emd'
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a, Applicative m)
    => (SiftResult v (n + 1) a -> m r)
    -> SiftCondition a
    -> SVG.Vector v (n + 1) a
    -> m (EMD v (n + 1) a)
emd'  cb sc = go id
  where
    go !imfs !v = cb res *> case res of
        SRResidual r -> pure $ EMD (imfs []) r
        SRIMF v' _   -> go (imfs . (v':)) (v - v')
      where
        res = sift sc v

data SiftResult v n a = SRResidual !(SVG.Vector v n a)
                      | SRIMF      !(SVG.Vector v n a) !Int   -- number of iterations

-- | Iterated sifting process.
sift
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SiftCondition a
    -> SVG.Vector v (n + 1) a
    -> SiftResult v (n + 1) a
sift sc = go 0
  where
    go !i !v = case sift' v of
      Nothing -> SRResidual v
      Just !v'
        | testCondition sc i v v' -> SRIMF v' i
        | otherwise               -> go (i + 1) v'

-- | Single sift
sift'
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SVG.Vector v (n + 1) a
    -> Maybe (SVG.Vector v (n + 1) a)
sift' v = go <$> envelopes v
  where
    go (mins, maxs) = SVG.zipWith3 (\x mi ma -> x - (mi + ma)/2) v mins maxs

-- | Returns cubic splines of local minimums and maximums.  Returns
-- 'Nothing' if there are not enough local minimum or maximums to create
-- the splines.
--
-- We add endpoints as both minima and maxima, to match behavior of
-- <http://www.mit.edu/~gari/CODE/HRV/emd.m>
envelopes
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SVG.Vector v (n + 1) a
    -> Maybe (SVG.Vector v (n + 1) a, SVG.Vector v (n + 1) a)
-- envelopes v = traceShow (M.size mins, M.size maxs) $ do
envelopes v = do
    guard $ (M.size mins + M.size maxs) >= 5
    (,) <$> splineAgainst mins
        <*> splineAgainst maxs
  where
    (mins, maxs) = extrema v

-- | Returns local minimums and maximums.
--
-- We add endpoints as both minima and maxima, to match behavior of
-- <http://www.mit.edu/~gari/CODE/HRV/emd.m>
extrema
    :: (VG.Vector v a, Fractional a, Ord a)
    => SVG.Vector v n a
    -> (M.Map (Finite n) a, M.Map (Finite n) a)
extrema = uncurry process
        . SVG.ifoldl' (uncurry go) (Nothing, ([], []))
  where
    go = \case
      Nothing            -> \_ i x -> (Just (x, EQ, i), ([(i,x)],[(i,x)]))
      -- Nothing            -> \_ i x -> (Just (x, EQ, i), ([],[]))
      Just (!y, !d, !i') -> \mms@(!mins, !maxs) !i !x ->
        let d'       = (x - y) `compare` 0
            newD     = case d' of
                         LT -> LT
                         EQ -> d
                         GT -> GT
            newLocal = (i', y)
            mms'     = case (d, d') of
              (LT, LT) -> mms
              (LT, EQ) -> (mins, maxs)
              (LT, GT) -> (newLocal : mins, maxs)
              (EQ, _ ) -> mms
              (GT, LT) -> (mins, newLocal : maxs)
              (GT, EQ) -> (mins, maxs)
              (GT, GT) -> mms
        in  (Just (x, newD, i), mms')
    process = \case
      -- Just (!y, _, !i') -> join bimap (M.fromList . ((i',y):))
      -- Nothing           -> join bimap M.fromList
      Just (!y, _, !i') -> join bimap (M.fromDescList . ((i',y):))
      Nothing           -> join bimap M.fromDescList

splineAgainst
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => M.Map (Finite n) a
    -> Maybe (SVG.Vector v n a)
splineAgainst = fmap go . makeSpline . M.mapKeysMonotonic fromIntegral
  where
    go spline = SVG.generate (sampleSpline spline . fromIntegral)
