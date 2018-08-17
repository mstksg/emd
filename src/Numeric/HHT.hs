{-# LANGUAGE BangPatterns                             #-}
{-# LANGUAGE DataKinds                                #-}
{-# LANGUAGE DeriveGeneric                            #-}
{-# LANGUAGE FlexibleContexts                         #-}
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
    HHT(..), HHTLine(..)
  , hhtEmd
  , hht
  , hhtSpectrum, hhtSpectrumT
  , marginal, instantaneousEnergy, degreeOfStationarity
  , expectedFreq
  -- ** Options
  , EMDOpts(..), defaultEO, BoundaryHandler(..), SiftCondition(..), defaultSC, SplineEnd(..)
  -- * Hilbert transforms (internal usage)
  , hilbert
  , hilbertIm
  , hilbertMagFreq
  ) where

import           Data.Complex
import           Data.Finite
import           Data.Fixed
import           Data.Foldable
import           Data.Proxy
import           Data.Semigroup
import           GHC.Generics              (Generic)
import           GHC.TypeNats
import           Numeric.EMD
import qualified Data.Binary               as Bi
import qualified Data.Map                  as M
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
  deriving (Show, Eq, Ord, Generic)

-- | @since 0.1.3.0
instance (VG.Vector v a, KnownNat n, Bi.Binary (v a)) => Bi.Binary (HHTLine v n a) where
    put HHTLine{..} = Bi.put (SVG.fromSized hlMags )
                   *> Bi.put (SVG.fromSized hlFreqs)
    get = do
      Just hlMags  <- SVG.toSized <$> Bi.get
      Just hlFreqs <- SVG.toSized <$> Bi.get
      pure HHTLine{..}

-- | A Hilbert-Huang Transform.  An @'HHT' v n a@ is a Hilbert-Huang
-- transform of an @n@-item time series of items of type @a@ represented
-- using vector @v@.
--
-- Create using 'hht' or 'hhtEmd'.
newtype HHT v n a = HHT { hhtLines :: [HHTLine v n a] }
  deriving (Show, Eq, Ord, Generic)

-- | @since 0.1.3.0
instance (VG.Vector v a, KnownNat n, Bi.Binary (v a)) => Bi.Binary (HHT v n a)

-- | Directly compute the Hilbert-Huang transform of a given time series.
-- Essentially is a composition of 'hhtEmd' and 'emd'.  See 'hhtEmd' for
-- a more flexible version.
hht :: forall v n a. (VG.Vector v a, KnownNat n, RealFloat a)
    => EMDOpts a
    -> SVG.Vector v (n + 1) a
    -> HHT v n a
hht eo = hhtEmd . emd eo

-- | Compute the Hilbert-Huang transform from a given Empirical Mode
-- Decomposition.
hhtEmd
    :: forall v n a. (VG.Vector v a, KnownNat n, RealFloat a)
    => EMD v (n + 1) a
    -> HHT v n a
hhtEmd EMD{..} = HHT $ map go emdIMFs
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
--
-- See 'httSpectrumT' for a verson with a "transposed" result.
hhtSpectrum
    :: forall v n a k. (VG.Vector v a, KnownNat n, Ord k, Num a)
    => (a -> k)     -- ^ binning function.  takes rev/tick freq between 0 and 1.
    -> HHT v n a
    -> SV.Vector n (M.Map k a)
hhtSpectrum f = foldl' ((SV.zipWith . M.unionWith) (+)) (pure mempty) . map go . hhtLines
  where
    go :: HHTLine v n a -> SV.Vector n (M.Map k a)
    go HHTLine{..} = SV.generate $ \i ->
      M.singleton (f $ hlFreqs `SVG.index` i) (hlMags `SVG.index` i)

-- | A "transposed" vesion of 'hhtSpectrum'.  Compute the full
-- Hilbert-Huang Transform spectrum.  Returns a map of the/sparse/ power
-- time series associated with each frequency.  Each frequency contains
-- a @'M.Map' ('Finite' n) a@, a map assocating times (the @'Finite' n@) to
-- the power at that time at that frequency.
--
-- Takes a "binning" function to allow you to specify how specific you want
-- your frequencies to be.
--
-- Note its relationship to the marginal power spectrum, 'marginal':
--
-- @
-- 'marginal' f = 'fmap' 'sum' . 'hhtSpectrumT' f
-- @
--
-- Where 'marginal' acts like a fourier transform (giving the total power
-- over the entire time series for every frequency), 'hhtSpectrumT' can be
-- thought of as giving the /temporal/ localizations of all contributions
-- of power from a given frequency.
--
-- @since 0.1.4.0
hhtSpectrumT
    :: forall v n a k. (VG.Vector v a, KnownNat n, Ord k, Num a)
    => (a -> k)     -- ^ binning function.  takes rev/tick freq between 0 and 1.
    -> HHT v n a
    -> M.Map k (M.Map (Finite n) a)
hhtSpectrumT f = M.unionsWith (M.unionWith (+)) . concatMap go . hhtLines
  where
    go :: HHTLine v n a -> [M.Map k (M.Map (Finite n) a)]
    go HHTLine{..} = flip fmap (finites @n) $ \i ->
      M.singleton (f $ hlFreqs `SVG.index` i) $
        M.singleton i (hlMags `SVG.index` i)

-- | Compute the marginal spectrum given a Hilbert-Huang Transform. It is
-- similar to a Fourier Transform; it provides the "total power" over the
-- entire time series for each frequency component.
--
-- A binning function is accepted to allow you to specify how specific you
-- want your frequencies to be.
--
-- Note the relationship to 'hhtSpectrumT':
--
-- @
-- 'marginal' f = 'fmap' 'sum' . 'hhtSpectrumT' f
-- @
--
-- It sums over the entire time series for every frequency.
marginal
    :: forall v n a k. (VG.Vector v a, KnownNat n, Ord k, Num a)
    => (a -> k)     -- ^ binning function.  takes rev/tick freq between 0 and 1.
    -> HHT v n a
    -> M.Map k a
marginal f = M.unionsWith (+) . concatMap go . hhtLines
  where
    go :: HHTLine v n a -> [M.Map k a]
    go HHTLine{..} = flip fmap (finites @n) $ \i ->
      M.singleton (f $ hlFreqs `SVG.index` i) (hlMags `SVG.index` i)

-- | Returns the "expected value" of frequency at each time step,
-- calculated as a weighted average of all contributions at every frequency
-- at that time step.
--
-- Essentially selects the dominant frequency at each time step.
--
-- @since 0.1.4.0
expectedFreq
    :: forall v n a. (VG.Vector v a, KnownNat n, Fractional a)
    => HHT v n a
    -> SVG.Vector v n a
expectedFreq HHT{..} = SVG.generate $ \i -> weightedAverage . map (go i) $ hhtLines
  where
    go :: Finite n -> HHTLine v n a -> (a, a)
    go i HHTLine{..} = (hlFreqs `SVG.index` i, hlMags `SVG.index` i)

weightedAverage
    :: (Foldable t, Fractional a)
    => t (a, a)
    -> a
weightedAverage = uncurry (/) . foldl' go (0, 0)
  where
    go (!sx, !sw) (!x, !w) = (sx + x, sw + w)

-- | Compute the instantaneous energy of the time series at every step via
-- the Hilbert-Huang Transform.
instantaneousEnergy
    :: forall v n a. (VG.Vector v a, KnownNat n, Num a)
    => HHT v n a
    -> SVG.Vector v n a
instantaneousEnergy = sum . map (SVG.map (^ (2 :: Int)) . hlMags) . hhtLines

-- | Degree of stationarity, as a function of frequency.
degreeOfStationarity
    :: forall v n a k. (VG.Vector v a, KnownNat n, Ord k, Fractional a)
    => (a -> k)     -- ^ binning function.  takes rev/tick freq between 0 and 1.
    -> HHT v n a
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
