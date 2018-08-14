{-# LANGUAGE DataKinds                                #-}
{-# LANGUAGE ExistentialQuantification                #-}
{-# LANGUAGE FlexibleContexts                         #-}
{-# LANGUAGE RankNTypes                               #-}
{-# LANGUAGE RecordWildCards                          #-}
{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE TypeApplications                         #-}
{-# LANGUAGE TypeOperators                            #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise       #-}

-- |
-- Module      : Numeric.HHT
-- Copyright   : (c) Justin Le 2018
-- License     : BSD3
--
-- Maintainer  : justin@jle.im
-- Stability   : experimental
-- Portability : non-portable
--
-- Hilbert-Huang transform in pure Haskell.
--
-- The main data type is 'HHT', which can be generated using 'hht' or
-- 'hhtEmd'.  See "Numeric.EMD" for information on why this module uses
-- "sized vectors", and how to convert unsized vectors to sized vectors.
--
-- @since 0.1.2.0

module Numeric.HHT (
  -- * Hilbert-Huang Transform
    hhtEmd
  , hht
  , hhtSpectrum
  -- ** Computing properties
  , marginal, instantaneousEnergy, degreeOfStationarity
  , HHT(..), HHTLine(..)
  , EMDOpts(..), defaultEO, BoundaryHandler(..), SiftCondition(..), defaultSC, SplineEnd(..)
  -- ** SomeHHT
  , SomeHHT(..), someHht
  -- * Hilbert transforms (internal usage)
  , hilbert
  , hilbertIm
  , hilbertMagFreq
  ) where

import           Data.Complex
import           Data.Finite
import           Data.Fixed
import           Data.List
import           Data.Proxy
import           Data.Semigroup
import           GHC.TypeNats
import           Numeric.EMD
import qualified Data.Map                  as M
import qualified Data.Vector               as V
import qualified Data.Vector.Generic       as VG
import qualified Data.Vector.Generic.Sized as SVG
import qualified Data.Vector.Sized         as SV

-- | A Hilbert Trasnform of a given IMF, given as a "skeleton line".
data HHTLine v n a = HHTLine
    { -- | IMF HHT Magnitude as a time series
      hlMags  :: !(SVG.Vector v n a)
      -- | IMF HHT instantaneous frequency as a time series (between 0 and 1)
    , hlFreqs :: !(SVG.Vector v n a)
    }
  deriving (Show, Eq, Ord)

-- | A Hilbert-Huang Transform.  An @'HHT' v n i a@ is a Hilbert-Huang
-- transform of an @n@-item time series with @i@ intrinsic mode functions
-- (IMFs) of items of type @a@ represented using vector @v@.
--
-- Create using 'hht' or 'hhtEmd'.
newtype HHT v n i a = HHT { hhtLines :: SV.Vector i (HHTLine v n a) }

-- | A @'SomeHHT' v n a@ is a convenient wrapper for an @'HHT' v n i a@
-- when you don't care about @i@, the number of intrinsic mode functions.
--
-- Pattern match to reveal the 'HHT' for functions that expect an 'EMD';
-- however, anything you do with the 'HHT' must be valid for /all/
-- potential IMF counts.
data SomeHHT v n a = forall i. KnownNat i => SomeHHT { getSomeHHT :: HHT v n i a }

-- | Directly compute the Hilbert-Huang transform of a given time series.
-- Essentially is a composition of 'hhtEmd' and 'emd'.  See 'hhtEmd' for
-- a more flexible version.
--
-- Returns an existentially quanitified number of IMF functions inside
-- a callback.  See 'someHht' as a version that returns a wrapped data type
-- instead.
hht :: forall v n a r. (VG.Vector v a, KnownNat n, RealFloat a)
    => EMDOpts a
    -> SVG.Vector v (n + 1) a
    -> (forall i. KnownNat i => HHT v n i a -> r)
    -> r
hht eo v f = emd eo v (f . hhtEmd)

-- | A version of 'hht' that returns a 'SomeHHT'.
someHht
    :: forall v n a. (VG.Vector v a, KnownNat n, RealFloat a)
    => EMDOpts a
    -> SVG.Vector v (n + 1) a
    -> SomeHHT v n a
someHht eo v = hht eo v SomeHHT

-- | Compute the Hilbert-Huang transform from a given Empirical Mode
-- Decomposition.
--
-- The returning 'HHT' will always have the same number of IMFs as the
-- given 'EMD', and this is enforced in the type.
hhtEmd
    :: forall v n a i. (VG.Vector v a, KnownNat n, RealFloat a)
    => EMD v (n + 1) i a
    -> HHT v n i a
hhtEmd EMD{..} = HHT $ SV.map go emdIMFs
  where
    go i = HHTLine (SVG.init m) f
      where
        (m, f) = hilbertMagFreq i

-- | Compute the full Hilbert-Huang Transform spectrum.  At each timestep
-- is a sparse map of frequency components and their respective magnitudes.
-- Frequencies not in the map are considered to be zero.
--
-- Takes a "binning" function to allow you to specify how specific you want
-- your frequencies to be.
hhtSpectrum
    :: forall n a k i. (KnownNat n, Ord k, Num a)
    => (a -> k)     -- ^ binning function.  takes rev/tick freq between 0 and 1.
    -> HHT V.Vector n i a
    -> SV.Vector n (M.Map k a)
hhtSpectrum f = foldl' ((SV.zipWith . M.unionWith) (+)) (pure mempty) . SV.map go . hhtLines
  where
    go :: HHTLine V.Vector n a -> SV.Vector n (M.Map k a)
    go HHTLine{..} = SV.generate $ \i ->
      M.singleton (f $ hlFreqs `SVG.index` i) (hlMags `SVG.index` i)

-- | Compute the marginal spectrum given a Hilbert-Huang Transform.
-- A binning function is accepted to allow you to specify how specific you
-- want your frequencies to be.
marginal
    :: forall v n a k i. (VG.Vector v a, KnownNat n, Ord k, Num a)
    => (a -> k)     -- ^ binning function.  takes rev/tick freq between 0 and 1.
    -> HHT v n i a
    -> M.Map k a
marginal f = M.unionsWith (+) . concatMap go . hhtLines
  where
    go :: HHTLine v n a -> [M.Map k a]
    go HHTLine{..} = flip fmap (finites @n) $ \i ->
      M.singleton (f $ hlFreqs `SVG.index` i) (hlMags `SVG.index` i)

-- | Compute the instantaneous energy of the time series at every step via
-- the Hilbert-Huang Transform.
instantaneousEnergy
    :: forall v n a i. (VG.Vector v a, KnownNat n, Num a)
    => HHT v n i a
    -> SVG.Vector v n a
instantaneousEnergy = sum . SV.map (SVG.map (^ (2 :: Int)) . hlMags) . hhtLines

-- | Degree of stationarity, as a function of frequency.
degreeOfStationarity
    :: forall v n a k i. (VG.Vector v a, KnownNat n, Ord k, Fractional a)
    => (a -> k)     -- ^ binning function.  takes rev/tick freq between 0 and 1.
    -> HHT v n i a
    -> M.Map k a
degreeOfStationarity f h = M.unionsWith (+)
                         . concatMap go
                         . hhtLines
                         $ h
  where
    meanMarg = (/ fromIntegral (natVal (Proxy @n))) <$> marginal f h
    go :: HHTLine v n a -> [M.Map k a]
    go HHTLine{..} = flip fmap (finites @n) $ \i ->
        let fr = f $ hlFreqs `SVG.index` i
        in M.singleton fr $
              (1 - (hlMags `SVG.index` i / meanMarg M.! fr)) ^ (2 :: Int)

-- | Given a time series, return a time series of the /magnitude/ of the
-- hilbert transform and the /frequency/ of the hilbert transform, in units
-- of revolutions per tick.  Is only expected to taken in proper/legal
-- IMFs.
--
-- The frequency will always be between 0 and 1, since we can't determine
-- anything faster given the discretization, and we exclude negative values
-- as physically unmeaningful for an IMF.
hilbertMagFreq
    :: forall v n a. (VG.Vector v a, KnownNat n, RealFloat a)
    => SVG.Vector v (n + 1) a
    -> (SVG.Vector v (n + 1) a, SVG.Vector v n a)
hilbertMagFreq v = (hilbertMag, hilbertFreq)
  where
    v' = hilbertIm v
    hilbertMag   = SVG.zipWith (\x x' -> magnitude (x :+ x')) v v'
    hilbertPhase = SVG.zipWith (\x x' -> phase (x :+ x')) v v'
    hilbertFreq  = SVG.map ((`mod'` 1) . (/ (2 * pi))) $ SVG.tail hilbertPhase - SVG.init hilbertPhase

-- | Real part is original series and imaginary part is hilbert transformed
-- series.  Creates a "helical" form of the original series that rotates
-- along the complex plane.
--
-- Numerically assumes that the signal is zero everywhere outside of the
-- vector, instead of the periodic assumption taken by matlab's version.
hilbert
    :: forall v n a. (VG.Vector v a, VG.Vector v (Complex a), KnownNat n, Floating a)
    => SVG.Vector v n a
    -> SVG.Vector v n (Complex a)
hilbert v = SVG.zipWith (:+) v (hilbertIm v)

-- | Hilbert transformed series.  Essentially the same series, but
-- phase-shifted 90 degrees.  Is so-named because it is the "imaginary
-- part" of the proper hilbert transform, 'hilbert'.
--
-- Numerically assumes that the signal is zero everywhere outside of the
-- vector, instead of the periodic assumption taken by matlab's version.
hilbertIm
    :: forall v n a. (VG.Vector v a, KnownNat n, Floating a)
    => SVG.Vector v n a
    -> SVG.Vector v n a
hilbertIm v = SVG.generate $ \i -> getSum . foldMap (Sum . go i) $ finites @n
  where
    -- NOTE: Can be made faster using an FFT and iFFT combo
    go :: Finite n -> Finite n -> a
    go i j
        | even k    = 0
        | otherwise = 2 * (v `SVG.index` j) / pi / fromIntegral k
      where
        k :: Int
        k = fromIntegral i - fromIntegral j
