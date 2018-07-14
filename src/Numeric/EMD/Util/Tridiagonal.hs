{-# LANGUAGE DataKinds                                #-}
{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE TypeApplications                         #-}
{-# LANGUAGE TypeFamilies                             #-}
{-# LANGUAGE TypeOperators                            #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise       #-}

module Numeric.EMD.Util.Tridiagonal (
    solveTridiagonal
  ) where

import           Control.Applicative.Backwards
import           Control.Monad.ST
import           Data.Finite
import           Data.Foldable
import           GHC.TypeNats
import qualified Data.Vector.Generic               as VG
import qualified Data.Vector.Generic.Mutable.Sized as SMVG
import qualified Data.Vector.Generic.Sized         as SVG

-- | https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
solveTridiagonal
    :: forall v n a. (VG.Vector v a, KnownNat n, Fractional a, 1 <= n)
    => SVG.Vector v n       a       -- ^ Bottom diagonal of M
    -> SVG.Vector v (n + 1) a       -- ^ Main diagonal of M
    -> SVG.Vector v n       a       -- ^ Upper diagonal of M
    -> SVG.Vector v (n + 1) a       -- ^ y
    -> SVG.Vector v (n + 1) a       -- ^ x such that M x = y
solveTridiagonal as bs cs ds = runST $ do
    mxs <- SVG.thaw ds'
    forwards . for_ (consecFinites @n) $ \(i0, i1) -> Backwards $ do
      x1 <- SMVG.read mxs i1
      let sbr = cs' `SVG.index` i0 * x1
      SMVG.modify mxs (subtract sbr) (weaken i0)
    SVG.freeze mxs
  where
    cs' :: SVG.Vector v n a
    cs' = runST $ do
      mcs <- SVG.thaw cs
      SMVG.modify mcs (/ SVG.head bs) minBound
      for_ (consecFinites @(n - 1)) $ \(i0, i1) -> do
        c0 <- SMVG.read mcs (weaken i0)
        let dvr = bs `SVG.index` weaken i1
                - as `SVG.index` weaken i0 * c0
        SMVG.modify mcs (/ dvr) i1
      SVG.freeze mcs
    ds' :: SVG.Vector v (n + 1) a
    ds' = runST $ do
      mds <- SVG.thaw ds
      SMVG.modify mds (/ SVG.head bs) minBound
      for_ (consecFinites @n) $ \(i0, i1) -> do
        let c0 = cs' `SVG.index` i0
        d0 <- SMVG.read mds (weaken i0)
        let sbr = as `SVG.index` i0 * d0
            dvr = bs `SVG.index` i1
                - as `SVG.index` i0 * c0
        SMVG.modify mds ((/ dvr) . subtract sbr) i1
      SVG.freeze mds

-- TODO: could be optimized unsafely
consecFinites :: KnownNat n => [(Finite n, Finite (n + 1))]
consecFinites = map (\i -> (i, shift i)) finites

    -- ds' :: SVG.Vector v (n + 1) a
    -- ds' = _
    -- xs  = _

    -- runST $ do
    -- cs' <- SVG.thaw cs
    -- _
