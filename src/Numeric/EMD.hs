{-# LANGUAGE BangPatterns                             #-}
{-# LANGUAGE DeriveGeneric                            #-}
{-# LANGUAGE GADTs                                    #-}
{-# LANGUAGE LambdaCase                               #-}
{-# LANGUAGE RankNTypes                               #-}
{-# LANGUAGE RecordWildCards                          #-}
{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE TypeApplications                         #-}
{-# LANGUAGE TypeInType                               #-}
{-# LANGUAGE TypeOperators                            #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise       #-}

-- |
-- Module      : Numeric.EMD
-- Copyright   : (c) Justin Le 2018
-- License     : BSD3
--
-- Maintainer  : justin@jle.im
-- Stability   : experimental
-- Portability : non-portable
--
-- Empirical Mode Decomposition in pure Haskell.
--
-- Main interface is 'emd', with 'defaultEO'.  A tracing version that
-- outputs a log to stdout is also available, as 'emdTrace'.  This can be
-- used to help track down a specific IMF that might be taking more time
-- than desired.
--
-- This package uses "sized vectors" as its main interface, to ensure:
--
-- 1.  The resulting 'EMD' contains IMFs that are all the same length as
--     the input vector
-- 2.  We provide a vector of size of at least one.
--
-- There are many functions to convert unsized vectors to sized vectors in
-- "Data.Vector.Sized" and associated modules, including 'toSized' (for
-- when you know the size at compile-time) and 'withSized' (for when you
-- don't).
--

module Numeric.EMD (
  -- * Empirical Mode Decomposition
    emd
  , emdTrace
  , emd'
  , iemd
  , EMD(..)
  -- ** Configuration
  , EMDOpts(..), defaultEO
  , BoundaryHandler(..)
  , SiftCondition(..), SiftProjection(..), defaultSC
  , scEnergyDiff
  , SplineEnd(..)
  -- * Internal
  , sift, SiftResult(..)
  , envelopes
  ) where

import           Control.DeepSeq
import           Control.Monad.IO.Class
import           Data.Functor.Identity
import           Data.List
import           GHC.Generics              (Generic)
import           GHC.TypeNats
import           Numeric.EMD.Internal.Sift
import           Text.Printf
import qualified Data.Binary               as Bi
import qualified Data.Vector.Generic       as VG
import qualified Data.Vector.Generic.Sized as SVG

-- | An @'EMD' v n a@ is an Empirical Mode Decomposition of a time series
-- with @n@ items of type @a@ stored in a vector @v@.
--
-- The component-wise sum of 'emdIMFs' and 'emdResidual' should yield
-- exactly the original series (see 'iemd').
data EMD v n a = EMD { emdIMFs     :: ![SVG.Vector v n a]
                     , emdResidual :: !(SVG.Vector v n a)
                     }
  deriving (Show, Generic, Eq, Ord)

-- | @since 0.1.5.0
instance NFData (v a) => NFData (EMD v n a)

-- | @since 0.1.3.0
instance (VG.Vector v a, KnownNat n, Bi.Binary (v a)) => Bi.Binary (EMD v n a) where
    put EMD{..} = Bi.put (SVG.fromSized <$> emdIMFs)
               *> Bi.put (SVG.fromSized emdResidual)
    get = do
      Just emdIMFs     <- traverse SVG.toSized <$> Bi.get
      Just emdResidual <- SVG.toSized <$> Bi.get
      pure EMD{..}

-- | EMD decomposition of a given time series with a given sifting stop
-- condition.
--
-- Takes a sized vector to ensure that:
--
-- 1.  The resulting 'EMD' contains IMFs that are all the same length as
--     the input vector
-- 2.  We provide a vector of size of at least one.
emd :: (VG.Vector v a, KnownNat n, Floating a, Ord a)
    => EMDOpts a
    -> SVG.Vector v (n + 1) a
    -> EMD v (n + 1) a
emd eo = runIdentity . emd' (const (pure ())) eo

-- | 'emd', but tracing results to stdout as IMFs are found.  Useful for
-- debugging to see how long each step is taking.
emdTrace
    :: (VG.Vector v a, KnownNat n, Floating a, Ord a, MonadIO m)
    => EMDOpts a
    -> SVG.Vector v (n + 1) a
    -> m (EMD v (n + 1) a)
emdTrace = emd' $ \case
    SRResidual _ -> liftIO $ putStrLn "Residual found."
    SRIMF _ i    -> liftIO $ printf "IMF found (%d sifts)\n" i

-- | 'emd' with a callback for each found IMF.
emd'
    :: (VG.Vector v a, KnownNat n, Floating a, Ord a, Applicative m)
    => (SiftResult v (n + 1) a -> m r)
    -> EMDOpts a
    -> SVG.Vector v (n + 1) a
    -> m (EMD v (n + 1) a)
emd' cb eo = go id
  where
    go !imfs !v = cb res *> case res of
        SRResidual r -> pure $ EMD (imfs []) r
        SRIMF v' _   -> go (imfs . (v':)) (v - v')
      where
        res = sift eo v

-- | Collapse an 'EMD' back into its original time series.  Should be
-- a left-inverse to 'emd': using 'iemd' on the result of 'emd' should give
-- back the original vector.
--
-- @since 0.1.5.0
iemd
    :: (VG.Vector v a, Num a)
    => EMD v n a
    -> SVG.Vector v n a
iemd EMD{..} = foldl' (SVG.zipWith (+)) emdResidual emdIMFs
