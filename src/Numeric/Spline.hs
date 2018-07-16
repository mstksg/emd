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

import           Data.Finite
import           Data.Proxy
import           Data.Type.Equality
import           GHC.TypeLits.Compare
import           GHC.TypeNats
import           Numeric.EMD.Util.Tridiagonal
import qualified Data.Map                         as M
import qualified Data.Vector.Sized                as SV

data SplineCoef a = SC { scAlpha  :: !a      -- ^ a
                       , scBeta   :: !a      -- ^ b
                       , scGamma0 :: !a      -- ^ y_{i-1}
                       , scGamma1 :: !a      -- ^ y_i
                       , scDelta  :: !a      -- ^ x_i - x_{i-1}
                       }
  deriving Show

data Spline a = Spline { splineHead :: !(a, SplineCoef a)
                       , splineTail :: !(M.Map a (SplineCoef a))
                       }
  deriving Show

runSplineCoef :: Fractional a => a -> SplineCoef a -> a -> a
runSplineCoef x0 (SC α β γ0 γ1 δ) x = q * γ0
                                    + t * γ1
                                    + t * q * (q * α + t * β)
  where
    t = (x - x0) / δ
    q = 1 - t

sampleSpline
    :: (Fractional a, Ord a)
    => Spline a
    -> a
    -> a
sampleSpline Spline{..} x = case x `M.lookupLE` splineTail of
    Nothing -> case splineHead of
      (x0, sc) -> runSplineCoef x0 sc x
    Just (x0, sc) -> runSplineCoef x0 sc x

-- | <https://en.wikipedia.org/wiki/Spline_interpolation#Interpolation_using_natural_cubic_spline>
makeSpline
    :: forall a. (Ord a, Fractional a)
    => M.Map a a
    -> Maybe (Spline a)
makeSpline ps = do
    (xy0, ps') <- M.minViewWithKey ps
    SV.withSizedList (M.toList ps') $ \(xsys :: SV.Vector n (a, a)) -> do
      Refl <- Proxy @1 `isLE` Proxy @n
      let xs, ys :: SV.Vector (n + 1) a
          (xs, ys) = SV.unzip $ xy0 `SV.cons` xsys
          dxs, dys :: SV.Vector n a
          dxs = SV.tail xs - SV.init xs
          rdxs :: SV.Vector n a
          rdxs = recip dxs
          rdxssq :: SV.Vector n a
          rdxssq = rdxs * rdxs
          dys  = SV.tail ys - SV.init ys
          dydxssq = dys * rdxssq
          mainDiag :: SV.Vector (n - 1) a
          mainDiag = SV.zipWith (\rdx0 rdx1 -> 2 * ( rdx0 + rdx1 ))
                        (SV.init rdxs)
                        (SV.tail rdxs)
          lowerDiag :: SV.Vector (n - 1) a
          lowerDiag = SV.take rdxs
          upperDiag :: SV.Vector (n - 1) a
          upperDiag = SV.tail rdxs
          rhs :: SV.Vector (n - 1) a
          rhs = SV.zipWith (\dydxsq0 dydxsq1 -> 3 * (dydxsq0 + dydxsq1))
                        (SV.init dydxssq)
                        (SV.tail dydxssq)
          -- TODO: allow specifying end conditions
          firstRow :: (a, a)     -- main, upper
          firstRow = ( rdxssq `SV.index` minBound + rdx12
                     , rdxssq `SV.index` minBound
                     + rdxssq `SV.index` shift minBound
                     + 2 * rdx12
                     )
            where
              rdx12  = rdxs `SV.index` minBound * rdxs `SV.index` shift minBound
          lastRow :: (a, a)         -- lower, main
          lastRow = ( - (rdxssq `SV.index` weaken maxBound)
                      - (rdxssq `SV.index` maxBound)
                      - 2 * rdx12
                    , - rdxssq `SV.index` maxBound - rdx12
                    )
            where
              rdx12  = rdxs `SV.index` maxBound * rdxs `SV.index` weaken maxBound
          endRhs :: (a, a)
          endRhs = ( 2 * (dydxssq `SV.index` minBound) * (rdxs `SV.index` minBound)
                   + 3 * (dydxssq `SV.index` minBound) * (rdxs `SV.index` shift minBound)
                   + (dydxssq `SV.index` shift minBound) * (rdxs `SV.index` shift minBound)
                   , - (dydxssq `SV.index` weaken maxBound) * (rdxs `SV.index` weaken maxBound)
                     - 3 * (dydxssq `SV.index` maxBound) * (rdxs `SV.index` weaken maxBound)
                     - 2 * (dydxssq `SV.index` maxBound) * (rdxs `SV.index` maxBound)
                   )
      solution <- solveTridiagonal (lowerDiag `SV.snoc` fst lastRow)
                                   (fst firstRow `SV.cons` mainDiag `SV.snoc` snd lastRow)
                                   (snd firstRow `SV.cons` upperDiag)
                                   (fst endRhs `SV.cons` rhs `SV.snoc` snd endRhs)
      let as :: SV.Vector n a
          as = SV.zipWith3 (\k dx dy -> k * dx - dy) (SV.init solution) dxs dys
          bs :: SV.Vector n a
          bs = SV.zipWith3 (\k dx dy -> - k * dx + dy) (SV.tail solution) dxs dys
          coefs :: SV.Vector n (a, SplineCoef a)
          coefs = SV.zipWith6 (\x α β γ0 γ1 δ -> (x, SC α β γ0 γ1 δ))
                    (SV.init xs) as bs (SV.init ys) (SV.tail ys) dxs

      pure $ Spline (SV.head coefs)
                    (M.fromAscList . SV.toList . SV.tail $ coefs)
