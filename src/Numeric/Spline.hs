{-# LANGUAGE ApplicativeDo                            #-}
{-# LANGUAGE DataKinds                                #-}
{-# LANGUAGE GADTs                                    #-}
{-# LANGUAGE RecordWildCards                          #-}
{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE TypeApplications                         #-}
{-# LANGUAGE TypeOperators                            #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise       #-}

module Numeric.Spline (
    Spline
  , makeSpline
  , sampleSpline
  ) where

import           Control.Applicative.Backwards
import           Control.Monad.ST
import           Data.Finite
import           Data.Foldable
import           Data.Proxy
import           Data.Type.Equality
import           GHC.TypeLits.Compare
import           GHC.TypeNats
import qualified Data.Map                      as M
import qualified Data.Vector.Mutable.Sized     as SMV
import qualified Data.Vector.Sized             as SV

data SplineCoef a = SC { scAlpha :: !a
                       , scBeta  :: !a
                       , scGamma :: !a
                       , scDelta :: !a
                       }
  deriving Show

data Spline a = Spline { splineHead :: !(a, SplineCoef a)
                       , splineTail :: !(M.Map a (SplineCoef a))
                       }
  deriving Show

runSplineCoef :: Num a => a -> SplineCoef a -> a -> a
runSplineCoef x0 (SC α β γ δ) x = α + sum (zipWith (*) [β, γ, δ] pows)
  where
    δx   = x - x0
    pows = iterate (* δx) δx

sampleSpline
    :: (Num a, Ord a)
    => Spline a
    -> a
    -> a
sampleSpline Spline{..} x = case x `M.lookupLE` splineTail of
    Nothing -> case splineHead of
      (x0, sc) -> runSplineCoef x0 sc x
    Just (x0, sc) -> runSplineCoef x0 sc x

-- | Computes a 'Spline' based on x-y coordinates, based on
-- <https://en.wikipedia.org/wiki/Spline_(mathematics)#Algorithm_for_computing_natural_cubic_splines>
--
-- Returns 'Nothing' if given an empty map.
makeSpline
    :: forall a. (Ord a, Fractional a)
    => M.Map a a
    -> Maybe (Spline a)
makeSpline ps = SV.withSizedList (M.toList ps) $ \(xsys :: SV.Vector n (a,a)) ->
    go xsys <$> Proxy @2 `isLE` Proxy @n
  where
    go  :: forall n. KnownNat n
        => SV.Vector n (a, a)
        -> ((2 <=? n) :~: 'True)
        -> Spline a
    go xsys Refl = runST $ do
        βs <- SMV.unsafeNew @(n - 1)
        γs <- SMV.unsafeNew @n
        δs <- SMV.unsafeNew @(n - 1)
        λs <- SMV.unsafeNew @n
        μs <- SMV.unsafeNew @(n - 1)
        ζs <- SMV.unsafeNew @n

        SMV.write λs 0 1
        SMV.write μs 0 0
        SMV.write ζs 0 0
        for_ (finites @(n - 2)) $ \i0 -> do
          let i1 = shift i0
              i2 = shift i1
              h0 = hs `SV.index` weaken i0
          μ0 <- SMV.read μs (weaken  i0)
          ζ0 <- SMV.read ζs (weakenN i0)
          let λ = 2 * ( xs `SV.index` i2
                      - xs `SV.index` weakenN i0
                      )
                - h0 * μ0
              μ = hs `SV.index` i1 / λ
              ζ = (as `SV.index` i0 - h0 * ζ0) / λ
          SMV.write λs (weaken i1) λ
          SMV.write μs i1          μ
          SMV.write ζs (weaken i1 :: Finite n) ζ
        SMV.write λs maxBound 1
        SMV.write γs maxBound 0
        SMV.write ζs maxBound 0

        forwards . for_ (finites @(n - 1)) $ \i1 -> Backwards $ do
          let i2 = shift i1
          γ2 <- SMV.read γs i2
          ζ  <- SMV.read ζs (weaken i1)
          μ  <- SMV.read μs i1
          let h = hs `SV.index` i1
              γ = ζ - μ * γ2
              β = (ys `SV.index` i2 - ys `SV.index` weaken i1) / h
                - h * (γ2 + 2 * γ) / 3
              δ = (γ2 - γ) / 3 / h
          SMV.write γs (weaken i1) γ
          SMV.write βs i1          β
          SMV.write δs i1          δ

        coefs <- SV.generateM @(n - 1) $ \i ->
          SC <$> pure (ys `SV.index` weaken i)
             <*> SMV.read βs i
             <*> SMV.read γs (weaken i)
             <*> SMV.read δs i

        let res :: SV.Vector (n - 1) (a, SplineCoef a)
            res = SV.zip (SV.init xs) coefs
        pure $ Spline (SV.head @(n - 2) res)
                      (M.fromAscList . SV.toList . SV.tail @(n - 2) $ res)

      where
        xs :: SV.Vector n a
        ys :: SV.Vector n a
        (xs, ys) = SV.unzip xsys
        hs :: SV.Vector (n - 1) a
        hs = SV.generate $ \i -> xs `SV.index` shift i
                               - xs `SV.index` weaken i
        as :: SV.Vector (n - 2) a
        as = SV.generate $ \i0 ->
          let i1 = shift i0
              i2 = shift i1
              y0 = ys `SV.index` weakenN i0
              y1 = ys `SV.index` weaken  i1
              y2 = ys `SV.index` i2
              h0 = hs `SV.index` weaken  i0
              h1 = hs `SV.index` i1
          in  3 * ((y2 - y1) / h1 - (y1 - y0) / h0)
