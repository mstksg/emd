{-# LANGUAGE BangPatterns                             #-}
{-# LANGUAGE GADTs                                    #-}
{-# LANGUAGE LambdaCase                               #-}
{-# LANGUAGE RecordWildCards                          #-}
{-# LANGUAGE ScopedTypeVariables                      #-}
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
-- Empirical Mode Decomposition (Hilbert-Huang Transform) in pure Haskell.
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
-- However, for convenience, "Numeric.EMD.Unsized" is provided with an
-- unsafe unsized interface.
--

module Numeric.EMD (
  -- * EMD (Hilbert-Huang Transform)
    emd
  , emdTrace
  , emd'
  , EMD(..)
  , EMDOpts(..), defaultEO, SiftCondition(..), defaultSC, SplineEnd(..)
  -- * Internal
  , sift, SiftResult(..)
  , envelopes
  ) where

import           Control.Monad.IO.Class
import           Data.Finite
import           Data.Functor.Identity
import           GHC.TypeNats
import           Numeric.EMD.Internal.Extrema
import           Numeric.EMD.Internal.Spline
import           Text.Printf
import qualified Data.Map                     as M
import qualified Data.Vector.Generic          as VG
import qualified Data.Vector.Generic.Sized    as SVG

-- | Options for EMD composition.
data EMDOpts a = EO { eoSiftCondition :: SiftCondition a  -- ^ stop condition for sifting
                    , eoSplineEnd     :: SplineEnd a      -- ^ end conditions for envelope splines
                    , eoClampEnvelope :: Bool             -- ^ if 'True', use time series endpoints as part of min and max envelopes
                    }
  deriving (Show, Eq, Ord)

-- | Default 'EMDOpts'
defaultEO :: Fractional a => EMDOpts a
defaultEO = EO { eoSiftCondition = defaultSC
               , eoSplineEnd     = SENotAKnot
               , eoClampEnvelope = True
               }


-- | Stop conditions for sifting process
data SiftCondition a = SCStdDev a         -- ^ Stop using standard "SD" method
                     | SCTimes Int        -- ^ Stop after a fixed number of iterations
                     | SCOr (SiftCondition a) (SiftCondition a)   -- ^ one or the other
                     | SCAnd (SiftCondition a) (SiftCondition a)  -- ^ both conditions met
  deriving (Show, Eq, Ord)

-- | Default 'SiftCondition'
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

-- | An @'EMD' v n a@ is a Hilbert-Huang transform of a time series with
-- @n@ items of type @a@ stored in a vector @v@.
data EMD v n a = EMD { emdIMFs     :: ![SVG.Vector v n a]
                     , emdResidual :: !(SVG.Vector v n a)
                     }
  deriving Show

-- | EMD decomposition (Hilbert-Huang Transform) of a given time series
-- with a given sifting stop condition.
--
-- Takes a sized vector to ensure that:
--
-- 1.  The resulting 'EMD' contains IMFs that are all the same length as
--     the input vector
-- 2.  We provide a vector of size of at least one.
emd :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => EMDOpts a
    -> SVG.Vector v (n + 1) a
    -> EMD v (n + 1) a
emd eo = runIdentity . emd' (const (pure ())) eo

-- | 'emd', but tracing results to stdout as IMFs are found.  Useful for
-- debugging to see how long each step is taking.
emdTrace
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a, MonadIO m)
    => EMDOpts a
    -> SVG.Vector v (n + 1) a
    -> m (EMD v (n + 1) a)
emdTrace = emd' $ \case
    SRResidual _ -> liftIO $ putStrLn "Residual found."
    SRIMF _ i    -> liftIO $ printf "IMF found (%d iterations)\n" i

-- | 'emd' with a callback for each found IMF.
emd'
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a, Applicative m)
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

-- | The result of a sifting operation.  Each sift either yields
-- a residual, or a new IMF.
data SiftResult v n a = SRResidual !(SVG.Vector v n a)
                      | SRIMF      !(SVG.Vector v n a) !Int   -- ^ number of iterations

-- | Iterated sifting process, used to produce either an IMF or a residual.
sift
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => EMDOpts a
    -> SVG.Vector v (n + 1) a
    -> SiftResult v (n + 1) a
sift EO{..} = go 1
  where
    go !i !v = case sift' eoSplineEnd eoClampEnvelope v of
      Nothing -> SRResidual v
      Just !v'
        | testCondition eoSiftCondition i v v' -> SRIMF v' i
        | otherwise                            -> go (i + 1) v'

-- | Single sift
sift'
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SplineEnd a
    -> Bool
    -> SVG.Vector v (n + 1) a
    -> Maybe (SVG.Vector v (n + 1) a)
sift' se cl v = go <$> envelopes se cl v
  where
    go (mins, maxs) = SVG.zipWith3 (\x mi ma -> x - (mi + ma)/2) v mins maxs

-- | Returns cubic splines of local minimums and maximums.  Returns
-- 'Nothing' if there are not enough local minimum or maximums to create
-- the splines.
envelopes
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SplineEnd a
    -> Bool
    -> SVG.Vector v (n + 1) a
    -> Maybe (SVG.Vector v (n + 1) a, SVG.Vector v (n + 1) a)
envelopes se cl xs = (,) <$> splineAgainst se mins'
                         <*> splineAgainst se maxs'
  where
    minMax = M.fromList [(minBound, SVG.head xs), (maxBound, SVG.last xs)]
    (mins,maxs) = extrema xs
    (mins', maxs')
      | cl        = (mins `M.union` minMax, maxs `M.union` minMax)
      | otherwise = (mins, maxs)

-- | Build a splined vector against a map of control points.
splineAgainst
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SplineEnd a
    -> M.Map (Finite n) a
    -> Maybe (SVG.Vector v n a)
splineAgainst se = fmap go . makeSpline se . M.mapKeysMonotonic fromIntegral
  where
    go spline = SVG.generate (sampleSpline spline . fromIntegral)
