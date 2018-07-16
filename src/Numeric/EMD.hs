{-# LANGUAGE BangPatterns                             #-}
{-# LANGUAGE GADTs                                    #-}
{-# LANGUAGE LambdaCase                               #-}
{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE TypeApplications                         #-}
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

import           Control.Applicative
import           Control.Monad
import           Control.Monad.IO.Class
import           Data.Bifunctor
import           Data.Bitraversable
import           Data.Finite
import           Data.Functor.Identity
import           Data.List
import           GHC.TypeNats
import           Numeric.Spline
import           Text.Printf
import qualified Data.Map                  as M
import qualified Data.Set                  as S
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
      let sd = SVG.sum $ SVG.zipWith (\x x' -> (x-x')^(2::Int) / (x^(2::Int) + eps)) v v'
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
    -> SVG.Vector v (n + 2) a
    -> EMD v (n + 2) a
emd sc = runIdentity . emd' (const (pure ())) sc

-- | 'emd', but tracing results as IMFs are found.
emdTrace
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a, MonadIO m)
    => SiftCondition a
    -> SVG.Vector v (n + 2) a
    -> m (EMD v (n + 2) a)
emdTrace = emd' $ \case
    SRResidual _ -> liftIO $ putStrLn "Residual found!"
    SRIMF _ i    -> liftIO $ printf "IMF found (%d iterations)\n" i

-- | 'emd' with an optional callback for each found IMF.
emd'
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a, Applicative m)
    => (SiftResult v (n + 2) a -> m r)
    -> SiftCondition a
    -> SVG.Vector v (n + 2) a
    -> m (EMD v (n + 2) a)
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
    -> SVG.Vector v (n + 2) a
    -> SiftResult v (n + 2) a
sift sc = go 1
  where
    go !i !v = case sift' v of
      Nothing -> SRResidual v
      Just !v'
        | testCondition sc i v v' -> SRIMF v' i
        | otherwise               -> go (i + 1) v'

-- | Single sift
sift'
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SVG.Vector v (n + 2) a
    -> Maybe (SVG.Vector v (n + 2) a)
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
    => SVG.Vector v (n + 2) a
    -> Maybe (SVG.Vector v (n + 2) a, SVG.Vector v (n + 2) a)
envelopes = join bitraverse splineAgainst <=< extrema

-- | Returns local minimums and maximums.
--
-- We add endpoints as both minima and maxima, to match behavior of
-- <http://www.mit.edu/~gari/CODE/HRV/emd.m>
extrema
    :: forall v n a. (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SVG.Vector v (n + 2) a
    -> Maybe (M.Map (Finite (n + 2)) a, M.Map (Finite (n + 2)) a)
extrema xs = case altList (reverse optima) of
    alo@(ys@(y:_), zs@(z:_))
      | y > z     -> Just $ join bimap makeMap (zs, ys)
      | otherwise -> Just $ join bimap makeMap alo
    _             -> Nothing
  where
    dxs :: SVG.Vector v (n + 1) a
    dxs = SVG.tail xs - SVG.init xs
    optima :: [Finite (n + 1)]
    optima = foldl' go [] (finites @n)
      where
        go os i = case (dx0, dx1) of
            (EQ, _ ) -> weaken i : os
            (LT, GT) -> shift  i : os
            (GT, LT) -> shift  i : os
            (_ , _ ) -> os
          where
            dx0 = compare (dxs `SVG.index` weaken i) 0
            dx1 = compare (dxs `SVG.index` shift i) 0
    minMax :: S.Set (Finite (n + 2))
    minMax = S.fromAscList [minBound, maxBound]
    makeMap :: [Finite (n + 1)] -> M.Map (Finite (n + 2)) a
    makeMap = M.fromSet (xs `SVG.index`)
            . S.union minMax
            . S.fromAscList
            . map weaken

altList :: [a] -> ([a], [a])
altList []       = ([] , [])
altList [x]      = ([x], [])
altList (x:y:zs) = bimap (x:) (y:) . altList $ zs

splineAgainst
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => M.Map (Finite n) a
    -> Maybe (SVG.Vector v n a)
splineAgainst = fmap go . makeSpline . M.mapKeysMonotonic fromIntegral
  where
    go spline = SVG.generate (sampleSpline spline . fromIntegral)
