{-# LANGUAGE BangPatterns                             #-}
{-# LANGUAGE GADTs                                    #-}
{-# LANGUAGE LambdaCase                               #-}
{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE TypeInType                               #-}
{-# LANGUAGE TypeOperators                            #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise       #-}

module Numeric.EMD (
    emd
  , EMD(..)
  , SiftCondition(..), defaultSC
  , sift
  -- * Debug
  , extrema
  , envelopes
  ) where

-- import           Data.Kind
import           Control.Monad
import           Data.Bifunctor
import           Data.Bitraversable
import           Data.Finite
import           Data.Maybe
import           Data.Monoid
import           Debug.Trace
import           GHC.TypeNats
import           Numeric.Spline
import qualified Data.Map                  as M
import qualified Data.Vector.Generic       as VG
import qualified Data.Vector.Generic.Sized as SVG

-- data SiftCondition :: Type -> Type where
--     SCStdDev :: a -> SiftCondition
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

emd :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SiftCondition a
    -> SVG.Vector v (n + 1) a
    -> EMD v (n + 1) a
emd sc = go mempty
  where
    go !imfs !v = case sift sc v of
      SRResidual r -> trace "found Residual" $ EMD (appEndo imfs []  ) r
      SRIMF v'     -> trace "found IMF" $ go (imfs <> Endo (v':)) (v - v')

data SiftResult v n a = SRResidual !(SVG.Vector v n a)
                      | SRIMF      !(SVG.Vector v n a)

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
        | testCondition sc i v v' -> SRIMF v'
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
    guard $ (M.size mins + M.size maxs) > 4
    (,) <$> splineAgainst mins
        <*> splineAgainst maxs
  where
    (mins, maxs) = extrema v

-- | Returns local minimums and maximums.
--
-- We add endpoints as both minima and maxima, to match behavior of
-- <http://www.mit.edu/~gari/CODE/HRV/emd.m>
extrema
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SVG.Vector v n a
    -> (M.Map (Finite n) a, M.Map (Finite n) a)
extrema = uncurry process
        . SVG.ifoldl' (uncurry go) (Nothing, ([], []))
  where
    go = \case
      Nothing            -> \_ i x -> (Just (x, EQ, i), ([(i,x)],[(i,x)]))
      Just (!y, !d, !i') -> \mms@(!mins, !maxs) !i !x ->
        let d'       = (x - y) `compare` 0
            newLocal = (i', y)
            mms'     = case (d, d') of
              (LT, LT) -> mms
              (LT, EQ) -> (newLocal : mins, maxs)
              (LT, GT) -> (newLocal : mins, maxs)
              (EQ, _ ) -> mms
              (GT, LT) -> (mins, newLocal : maxs)
              (GT, EQ) -> (mins, newLocal : maxs)
              (GT, GT) -> mms
        in  (Just (x, d', i), mms')
    process = \case
      Just (!y, _, !i') -> join bimap (M.fromDescList . ((i',y):))
      Nothing           -> join bimap M.fromDescList
    -- process ()


splineAgainst
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => M.Map (Finite n) a
    -> Maybe (SVG.Vector v n a)
splineAgainst = fmap go . makeSpline . M.mapKeysMonotonic fromIntegral
  where
    go spline = SVG.generate (sampleSpline spline . fromIntegral)
