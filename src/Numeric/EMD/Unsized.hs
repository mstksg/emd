{-# LANGUAGE DataKinds                                #-}
{-# LANGUAGE GADTs                                    #-}
{-# LANGUAGE LambdaCase                               #-}
{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE TypeApplications                         #-}
{-# LANGUAGE TypeOperators                            #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise       #-}

-- |
-- Module      : Numeric.EMD.Unsized
-- Copyright   : (c) Justin Le 2018
-- License     : BSD3
--
-- Maintainer  : justin@jle.im
-- Stability   : experimental
-- Portability : non-portable
--
-- Interface of "Numeric.EMD" re-exported in a non-typesafe "unsized" form.
-- Can be more convenient in certain situations, but "Numeric.EMD" is
-- recommended and preferred.
--


module Numeric.EMD.Unsized (
  -- -- * EMD (Hilbert-Huang Transform)
    emd
  , emdTrace
  , emd'
  , EMD(..)
  , E.EMDOpts(..), E.defaultEO, E.SiftCondition(..), E.defaultSC, E.SplineEnd(..)
  -- -- * Internal
  , sift, SiftResult(..)
  -- , envelopes
  ) where

import           Control.Monad.IO.Class
import           Data.Proxy
import           Data.Type.Equality
import           GHC.TypeLits.Compare
import           GHC.TypeNats
import qualified Data.Vector.Generic       as VG
import qualified Data.Vector.Generic.Sized as SVG
import qualified Numeric.EMD               as E

-- | An @'EMD' v a@ is a Hilbert-Huang transform of a time series with
-- items of type @a@ stored in a vector @v@.
data EMD v a = EMD { emdIMFs     :: ![v a]
                   , emdResidual :: !(v a)
                   }
  deriving Show

-- | The result of a sifting operation.  Each sift either yields
-- a residual, or a new IMF.
data SiftResult v a = SRResidual !(v a)
                    | SRIMF      !(v a) !Int   -- ^ number of iterations


-- | EMD decomposition (Hilbert-Huang Transform) of a given time series
-- with a given sifting stop condition.
--
-- Returns 'Nothing' if given an empty vector.
--
-- See 'Numeric.EMD.emd' for a type-safe version with guaruntees on the
-- output vector sizes.
emd :: (VG.Vector v a, Fractional a, Ord a)
    => E.EMDOpts a
    -> v a
    -> Maybe (EMD v a)
emd eo v = SVG.withSized v $ \(v' :: SVG.Vector v n a) -> do
    Refl <- Proxy @1 `isLE` Proxy @n
    pure . convertEMD $ E.emd @_ @_ @(n - 1) eo v'

-- | 'emd', but tracing results to stdout as IMFs are found.  Useful for
-- debugging to see how long each step is taking.
--
-- Returns 'Nothing' if given an empty vector.
emdTrace
    :: (VG.Vector v a, Fractional a, Ord a, MonadIO m)
    => E.EMDOpts a
    -> v a
    -> m (Maybe (EMD v a))
emdTrace eo v = SVG.withSized v $ \(v' :: SVG.Vector v n a) ->
    case Proxy @1 `isLE` Proxy @n of
      Nothing   -> pure Nothing
      Just Refl -> Just . convertEMD <$> E.emdTrace @_ @_ @(n - 1) eo v'

-- | 'emd' with a callback for each found IMF.
--
-- Returns 'Nothing' if given an empty vector.
emd'
    :: (VG.Vector v a, Ord a, Fractional a, Applicative m)
    => (SiftResult v a -> m r)
    -> E.EMDOpts a
    -> v a
    -> m (Maybe (EMD v a))
emd' cb eo v = SVG.withSized v $ \(v' :: SVG.Vector v n a) ->
    case Proxy @1 `isLE` Proxy @n of
      Nothing   -> pure Nothing
      Just Refl -> Just . convertEMD <$> E.emd' @_ @_ @(n - 1) (cb . convertSR) eo v'

-- emd' cb eo = go id
--   where
--     go !imfs !v = cb res *> case res of
--         SRResidual r -> pure $ EMD (imfs []) r
--         SRIMF v' _   -> go (imfs . (v':)) (v - v')
--       where
--         res = sift eo v

-- | Iterated sifting process, used to produce either an IMF or a residual.
--
-- Returns 'Nothing' if given an empty vector.
sift
    :: (VG.Vector v a, Fractional a, Ord a)
    => E.EMDOpts a
    -> v a
    -> Maybe (SiftResult v a)
sift eo v = SVG.withSized v $ \(v' :: SVG.Vector v n a) -> do
    Refl <- Proxy @1 `isLE` Proxy @n
    pure $ convertSR . E.sift @_ @_ @(n - 1) eo $ v'

convertSR :: E.SiftResult v n a -> SiftResult v a
convertSR = \case
    E.SRResidual v -> SRResidual $ SVG.fromSized v
    E.SRIMF v i    -> SRIMF (SVG.fromSized v) i

convertEMD :: E.EMD v n a -> EMD v a
convertEMD (E.EMD is r) = EMD (SVG.fromSized <$> is) (SVG.fromSized r)

