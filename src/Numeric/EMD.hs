{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE LambdaCase   #-}
{-# LANGUAGE LambdaCase   #-}

module Numeric.EMD (
    sift
  , SiftCondition(..)
  ) where

import           Control.Monad.Trans.State
import           Control.Monad.Trans.Writer
import           Data.Bifunctor
import           Data.Finite
import           Data.Monoid
import           Data.Ord
import           GHC.TypeNats
import           Numeric.Spline
import qualified Data.Map                   as M
import qualified Data.Vector.Generic        as VG
import qualified Data.Vector.Generic.Sized  as SVG

data SiftCondition a = SCStdDev a

-- | 'True' if stop
testCondition
    :: (VG.Vector v a, Fractional a, Ord a)
    => SiftCondition a
    -> SVG.Vector v n a
    -> SVG.Vector v n a
    -> Bool
testCondition = \case
    SCStdDev t -> \v v' ->
      let sd = SVG.sum $ SVG.zipWith (\x x' -> (x-x')^(2::Int) / x^(2::Int)) v v'
      in  sd <= t

sift'
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SVG.Vector v n a
    -> SVG.Vector v n a
sift' v = v - ((mins + maxs) / 2)
  where
    (mins, maxs) = envelopes v

sift
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SiftCondition a
    -> SVG.Vector v n a
    -> SVG.Vector v n a
sift sc = go
  where
    go !v
        | testCondition sc v v' = v'
        | otherwise             = go v'
      where
        v' = sift' v

envelopes
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SVG.Vector v n a
    -> (SVG.Vector v n a, SVG.Vector v n a)
envelopes = bimap splineAgainst splineAgainst
          . snd
          . SVG.ifoldl' (uncurry go) (Nothing, (M.empty, M.empty))
  where
    go = \case
      Nothing            -> \mms i x ->  (Just (x, EQ, i), mms)
      Just (!y, !d, !i') -> \mms@(!mins, !maxs) !i !x ->
        let d'  = (x - y) `compare` 0
            mms' = case (d, d') of
              (LT, LT) -> mms
              (LT, EQ) -> (M.insert i' y mins, maxs)
              (LT, GT) -> (M.insert i' y mins, maxs)
              (EQ, _ ) -> mms
              (GT, LT) -> (mins, M.insert i' y maxs)
              (GT, EQ) -> (mins, M.insert i' y maxs)
              (GT, GT) -> mms
        in  (Just (x, d', i), mms')

splineAgainst
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => M.Map (Finite n) a
    -> SVG.Vector v n a
splineAgainst m = SVG.generate (sampleSpline spline . fromIntegral)
  where
    spline = makeSpline $ M.mapKeysMonotonic fromIntegral m

