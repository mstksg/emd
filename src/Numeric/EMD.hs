{-# LANGUAGE BangPatterns                             #-}
{-# LANGUAGE GADTs                                    #-}
{-# LANGUAGE LambdaCase                               #-}
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
  , EMD(..)
  , EMDOpts(..), defaultEO, BoundaryHandler(..), SiftCondition(..), defaultSC, SplineEnd(..)
  -- * Internal
  , sift, SiftResult(..)
  , envelopes
  ) where

import           Control.Monad
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
data EMDOpts a = EO { eoSiftCondition   :: SiftCondition a  -- ^ stop condition for sifting
                    , eoSplineEnd       :: SplineEnd a      -- ^ end conditions for envelope splines
                    , eoBoundaryHandler :: Maybe BoundaryHandler  -- ^ process for handling boundary
                    }
  deriving (Show, Eq, Ord)

data BoundaryHandler
    -- | Clamp envelope at end points (Matlab implementation)
    = BHClamp
    -- | Extend boundaries symmetrically
    | BHSymmetric
  deriving (Show, Eq, Ord)

    -- -- | Extend boundaries assuming global periodicity
    -- -- | BHPeriodic

-- | Default 'EMDOpts'
defaultEO :: Fractional a => EMDOpts a
defaultEO = EO { eoSiftCondition   = defaultSC
               , eoSplineEnd       = SENatural
               , eoBoundaryHandler = Just BHSymmetric
               }

-- | Stop conditions for sifting process
data SiftCondition a
    -- | Stop using standard SD method
    = SCStdDev !a
    -- | Stop after a fixed number of sifting iterations
    | SCTimes !Int
    -- | One or the other
    | SCOr (SiftCondition a) (SiftCondition a)
    -- | Stop when both conditions are met
    | SCAnd (SiftCondition a) (SiftCondition a)
  deriving (Show, Eq, Ord)

-- | Default 'SiftCondition'
defaultSC :: Fractional a => SiftCondition a
defaultSC = SCStdDev 0.3 `SCOr` SCTimes 100     -- R package uses SCTimes 20, Matlab uses no limit
-- defaultSC = SCStdDev 0.3

-- | 'True' if stop
testCondition
    :: (VG.Vector v a, Fractional a, Ord a)
    => SiftCondition a
    -> Int
    -> SVG.Vector v n a
    -> SVG.Vector v n a
    -> Bool
testCondition tc i v v' = go tc
  where
    sd = SVG.sum $ SVG.zipWith (\x x' -> (x-x')^(2::Int) / (x^(2::Int) + eps)) v v'
    go = \case
      SCStdDev t -> sd <= t
      SCTimes l  -> i >= l
      SCOr  f g  -> go f || go g
      SCAnd f g  -> go f && go g
    eps = 0.0000001

-- | An @'EMD' v n a@ is an Empirical Mode Decomposition of a time series
-- with @n@ items of type @a@ stored in a vector @v@.
--
-- The component-wise sum of 'emdIMFs' and 'emdResidual' should yield
-- exactly the original series.
data EMD v n a = EMD { emdIMFs     :: ![SVG.Vector v n a]
                     , emdResidual :: !(SVG.Vector v n a)
                     }
  deriving Show

-- | EMD decomposition of a given time series with a given sifting stop
-- condition.
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
    SRIMF _ i    -> liftIO $ printf "IMF found (%d sifts)\n" i

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
                      | SRIMF      !(SVG.Vector v n a) !Int   -- ^ number of sifting iterations

-- | Iterated sifting process, used to produce either an IMF or a residual.
sift
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => EMDOpts a
    -> SVG.Vector v (n + 1) a
    -> SiftResult v (n + 1) a
sift EO{..} = go 1
  where
    go !i !v = case sift' eoSplineEnd eoBoundaryHandler v of
      Nothing -> SRResidual v
      Just !v'
        | testCondition eoSiftCondition i v v' -> SRIMF v' i
        | otherwise                            -> go (i + 1) v'

-- | Single sift
sift'
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SplineEnd a
    -> Maybe BoundaryHandler
    -> SVG.Vector v (n + 1) a
    -> Maybe (SVG.Vector v (n + 1) a)
sift' se bh v = go <$> envelopes se bh v
  where
    go (mins, maxs) = SVG.zipWith3 (\x mi ma -> x - (mi + ma)/2) v mins maxs

-- | Returns cubic splines of local minimums and maximums.  Returns
-- 'Nothing' if there are not enough local minimum or maximums to create
-- the splines.
envelopes
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SplineEnd a
    -> Maybe BoundaryHandler
    -> SVG.Vector v (n + 1) a
    -> Maybe (SVG.Vector v (n + 1) a, SVG.Vector v (n + 1) a)
envelopes se bh xs = do
    when (bh == Just BHClamp) $ do
      guard (M.size mins > 1)
      guard (M.size maxs > 1)
    (,) <$> splineAgainst se emin mins
        <*> splineAgainst se emax maxs
  where
    -- minMax = M.fromList [(minBound, SVG.head xs), (maxBound, SVG.last xs)]
    (mins,maxs) = extrema xs
    (emin,emax) = case bh of
      Nothing  -> mempty
      Just bh' -> extendExtrema xs bh' (mins,maxs)
    --   | isJust bh = (mins `M.union` minMax, maxs `M.union` minMax)
    --   | otherwise = (mins, maxs)

extendExtrema
    :: forall v n a. (VG.Vector v a, KnownNat n)
    => SVG.Vector v (n + 1) a
    -> BoundaryHandler
    -> (M.Map (Finite (n + 1)) a, M.Map (Finite (n + 1)) a)
    -> (M.Map Int a, M.Map Int a)
    -- (M.Map (Finite (n + 1)) a, M.Map (Finite (n + 1)) a)
extendExtrema xs = \case
    BHClamp     -> const (firstLast, firstLast)
    BHSymmetric -> \(mins, maxs) ->
      let addFirst = case (flippedMin, flippedMax) of
              (Nothing      , Nothing      ) -> mempty
              -- first point is local maximum
              (Just (_,mn)  , Nothing      ) -> (mn        , firstPoint)
              -- first point is local minimum
              (Nothing      , Just (_,mx)  ) -> (firstPoint, mx        )
              (Just (mni,mn), Just (mxi,mx))
                | mni < mxi                  -> (mn        , firstPoint)
                | otherwise                  -> (firstPoint, mx        )
            where
              flippedMin = flip fmap (M.lookupMin mins) $ \(minIx, minVal) ->
                (minIx, M.singleton (negate (fromIntegral minIx)) minVal)
              flippedMax = flip fmap (M.lookupMin maxs) $ \(maxIx, maxVal) ->
                (maxIx, M.singleton (negate (fromIntegral maxIx)) maxVal)
          addLast = case (flippedMin, flippedMax) of
              (Nothing      , Nothing      ) -> mempty
              -- last point is local maximum
              (Just (_,mn)  , Nothing      ) -> (mn        , lastPoint )
              -- last point is local minimum
              (Nothing      , Just (_,mx)  ) -> (lastPoint , mx        )
              (Just (mni,mn), Just (mxi,mx))
                | mni > mxi                  -> (mn        , lastPoint )
                | otherwise                  -> (lastPoint , mx        )
            where
              flippedMin = flip fmap (M.lookupMax mins) $ \(minIx, minVal) ->
                (minIx, M.singleton (extendSym (fromIntegral minIx)) minVal)
              flippedMax = flip fmap (M.lookupMax maxs) $ \(maxIx, maxVal) ->
                (maxIx, M.singleton (extendSym (fromIntegral maxIx)) maxVal)
      in  addFirst `mappend` addLast
  where
    lastIx = fromIntegral $ maxBound @(Finite n)
    firstPoint = M.singleton 0 (SVG.head xs)
    lastPoint  = M.singleton lastIx (SVG.last xs)
    firstLast  = firstPoint `mappend` lastPoint
    extendSym i = 2 * lastIx - i

-- | Build a splined vector against a map of control points.
splineAgainst
    :: (VG.Vector v a, KnownNat n, Fractional a, Ord a)
    => SplineEnd a
    -> M.Map Int a              -- ^ extensions
    -> M.Map (Finite n) a
    -> Maybe (SVG.Vector v n a)
splineAgainst se ext = fmap go
                     . makeSpline se
                     . mappend (M.mapKeysMonotonic fromIntegral ext)
                     . M.mapKeysMonotonic fromIntegral
  where
    go spline = SVG.generate (sampleSpline spline . fromIntegral)
