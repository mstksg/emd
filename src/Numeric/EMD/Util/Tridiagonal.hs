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
import           Control.Monad
import           Control.Monad.ST
import           Control.Monad.Trans.Class
import           Control.Monad.Trans.Maybe
import           Data.Finite
import           Data.Foldable
import           GHC.TypeNats
import qualified Data.Vector.Generic               as VG
import qualified Data.Vector.Generic.Mutable.Sized as SMVG
import qualified Data.Vector.Generic.Sized         as SVG

-- | <https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm>
--
-- Will return 'Nothing' if the matrix is not invertible.  This will happen
-- if:
--
-- 1. The first item in the main diagonal is zero
-- 2. There is any i such that b_{i + 1} = a_i * c_i.  That is, an item in
-- the main diagonal is equal to the product of the off-diagonal elements
-- a row above it
-- 3. Another mystery condition!
solveTridiagonal
    :: forall v n a. (VG.Vector v a, KnownNat n, Fractional a, 1 <= n, Eq a)
    => SVG.Vector v n       a       -- ^ a: Bottom diagonal of M
    -> SVG.Vector v (n + 1) a       -- ^ b: Main diagonal of M
    -> SVG.Vector v n       a       -- ^ c: Upper diagonal of M
    -> SVG.Vector v (n + 1) a       -- ^ y
    -> Maybe (SVG.Vector v (n + 1) a) -- ^ x such that M x = y
solveTridiagonal as bs cs ds = runST $ runMaybeT $ do
    guard $ SVG.head bs /= 0
    cs' <- MaybeT . pure $ mcs'
    ds' <- MaybeT . pure $ mds'
    mxs <- lift $ SVG.thaw ds'
    forwards . for_ (consecFinites @n) $ \(i0, i1) -> Backwards $ do
      x1 <- lift $ SMVG.read mxs i1
      let sbr = cs' `SVG.index` i0 * x1
      lift $ SMVG.modify mxs (subtract sbr) (weaken i0)
    lift $ SVG.freeze mxs
  where
    mcs' :: Maybe (SVG.Vector v n a)
    mcs' = runST $ runMaybeT $ do
      mcs <- lift $ SVG.thaw cs
      lift $ SMVG.modify mcs (/ SVG.head bs) minBound
      for_ (consecFinites @(n - 1)) $ \(i0, i1) -> do
        c0 <- lift $ SMVG.read mcs (weaken i0)
        let dvr = bs `SVG.index` weaken i1
                - as `SVG.index` weaken i0 * c0
        guard $ dvr /= 0
        lift $ SMVG.modify mcs (/ dvr) i1
      lift $ SVG.freeze mcs
    mds' :: Maybe (SVG.Vector v (n + 1) a)
    mds' = runST $ runMaybeT $ do
      mds <- lift $ SVG.thaw ds
      cs' <- MaybeT . pure $ mcs'
      lift $ SMVG.modify mds (/ SVG.head bs) minBound
      for_ (consecFinites @n) $ \(i0, i1) -> do
        let c0 = cs' `SVG.index` i0
        d0 <- lift $ SMVG.read mds (weaken i0)
        let sbr = as `SVG.index` i0 * d0
            dvr = bs `SVG.index` i1
                - as `SVG.index` i0 * c0
        guard $ dvr /= 0
        lift $ SMVG.modify mds ((/ dvr) . subtract sbr) i1
      lift $ SVG.freeze mds

-- TODO: could be optimized unsafely
consecFinites :: KnownNat n => [(Finite n, Finite (n + 1))]
consecFinites = map (\i -> (i, shift i)) finites

    -- ds' :: SVG.Vector v (n + 1) a
    -- ds' = _
    -- xs  = _

    -- runST $ do
    -- cs' <- SVG.thaw cs
    -- _
