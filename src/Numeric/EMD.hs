{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE LambdaCase          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeInType          #-}

module Numeric.EMD (
    emd
  , EMD(..)
  , SiftCondition(..), defaultSC
  , sift
  ) where

-- import           Data.Kind
import           Control.Monad
import           Data.Bitraversable
import           Data.Finite
import           Data.Maybe
import           Data.Monoid
import           GHC.TypeNats
import           Numeric.Spline
import qualified Data.Map                  as M
import qualified Data.Vector.Generic       as VG
import qualified Data.Vector.Generic.Sized as SVG

-- data SiftCondition :: Type -> Type where
--     SCStdDev :: a -> SiftCondition
data SiftCondition a = SCStdDev a
                     | SCTimes Int
  deriving Show

defaultSC :: Fractional a => SiftCondition a
defaultSC = SCStdDev 0.3

-- | 'True' if stop
testCondition
    :: (VG.Vector v a, Fractional a, Ord a)
    => SiftCondition a
    -> SVG.Vector v n a
    -> SVG.Vector v n a
    -> Bool
testCondition = \case
    SCStdDev t -> \v v' ->
      let sd = SVG.sum $ SVG.zipWith (\x x' -> (x-x')^(2::Int) / (x^(2::Int) + eps)) v v'
      in  sd <= t
    SCTimes _  -> \_ _ -> True
  where
    eps = 0.0000001

data EMD v n a = EMD { emdIMFs     :: ![SVG.Vector v n a]
                     , emdResidual :: !(SVG.Vector v n a)
                     }
  deriving Show

emd :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SiftCondition a
    -> SVG.Vector v n a
    -> EMD v n a
emd sc = go mempty
  where
    go !imfs !v = case sift sc v of
      Nothing -> EMD (appEndo imfs []   ) v
      Just v' -> go  (imfs <> Endo (v':)) (v - v')

-- | Iterated sifting process
sift
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SiftCondition a
    -> SVG.Vector v n a
    -> Maybe (SVG.Vector v n a)
sift sc = go
  where
    go !v = nextStep <$> sift' v
      where
        nextStep v'
          | testCondition sc v v' = v'
          | otherwise             = fromMaybe v' (go v')

sift'
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SVG.Vector v n a
    -> Maybe (SVG.Vector v n a)
sift' v = go <$> envelopes v
  where
    go (mins, maxs) = SVG.zipWith3 (\x mi ma -> x - (mi + ma)/2) v mins maxs

-- | Returns cubic splines of local minimums and maximums.  Returns
-- 'Nothing' if there are not enough local minimum or maximums to create
-- the splines.
envelopes
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SVG.Vector v n a
    -> Maybe (SVG.Vector v n a, SVG.Vector v n a)
envelopes = join bitraverse (splineAgainst . M.fromDescList)
          . snd
          . SVG.ifoldl' (uncurry go) (Nothing, ([], []))
  where
    go = \case
      Nothing            -> \mms i x ->  (Just (x, EQ, i), mms)
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

splineAgainst
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => M.Map (Finite n) a
    -> Maybe (SVG.Vector v n a)
splineAgainst = fmap go . makeSpline . M.mapKeysMonotonic fromIntegral
  where
    go spline = SVG.generate (sampleSpline spline . fromIntegral)
