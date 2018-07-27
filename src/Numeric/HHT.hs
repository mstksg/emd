{-# LANGUAGE DataKinds                                #-}
{-# LANGUAGE FlexibleContexts                         #-}
{-# LANGUAGE RecordWildCards                          #-}
{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE TypeApplications                         #-}
{-# LANGUAGE TypeOperators                            #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise       #-}

module Numeric.HHT (
    hhtEmd
  , hht
  , HHT(..), HHTLine(..)
  , EMDOpts(..), defaultEO, SiftCondition(..), defaultSC, SplineEnd(..)
  -- * Hilbert transforms (internal usage)
  , hilbert
  , hilbert'
  ) where

import           Data.Complex
import           Data.Finite
import           Data.Fixed
import           Data.Semigroup
import           GHC.TypeNats
import           Numeric.EMD
import qualified Data.Vector.Generic       as VG
import qualified Data.Vector.Generic.Sized as SVG


data HHTLine v n a = HHTLine { hlMags  :: !(SVG.Vector v n a)
                             , hlFreqs :: !(SVG.Vector v n a)
                             }

newtype HHT v n a = HHT { hhtLines :: [HHTLine v n a] }

hht :: forall v n a. (VG.Vector v a, KnownNat n, RealFloat a)
    => EMDOpts a
    -> SVG.Vector v (n + 1) a
    -> HHT v n a
hht eo = hhtEmd . emd eo

hhtEmd
    :: forall v n a. (VG.Vector v a, KnownNat n, RealFloat a)
    => EMD v (n + 1) a
    -> HHT v n a
hhtEmd EMD{..} = HHT $ map go emdIMFs
  where
    go i = HHTLine (SVG.init m) f
      where
        (m, f) = hilbertMagFreq i

hilbertMagFreq
    :: forall v n a. (VG.Vector v a, KnownNat n, RealFloat a)
    => SVG.Vector v (n + 1) a
    -> (SVG.Vector v (n + 1) a, SVG.Vector v n a)
hilbertMagFreq v = (hilbertMag, hilbertFreq)
  where
    v' = hilbert' v
    hilbertMag   = SVG.zipWith (\x x' -> magnitude (x :+ x')) v v'
    hilbertPhase = SVG.zipWith (\x x' -> phase (x :+ x')) v v'
    hilbertFreq  = SVG.map wrap $ SVG.tail hilbertPhase - SVG.init hilbertPhase
    wrap          = subtract pi . (`mod'` (2 * pi)) . (+ pi)

-- | Real part is original series and imaginary part is hilbert transformed
-- series.
hilbert
    :: forall v n a. (VG.Vector v a, VG.Vector v (Complex a), KnownNat n, Floating a)
    => SVG.Vector v n a
    -> SVG.Vector v n (Complex a)
hilbert v = SVG.zipWith (:+) v (hilbert' v)

-- | Can be made faster using an FFT and iFFT combo
hilbert'
    :: forall v n a. (VG.Vector v a, KnownNat n, Floating a)
    => SVG.Vector v n a
    -> SVG.Vector v n a
hilbert' v = SVG.generate $ \i -> getSum . foldMap (Sum . go i) $ finites @n
  where
    go :: Finite n -> Finite n -> a
    go i j
        | even k    = 0
        | otherwise = 2 * (v `SVG.index` j) / pi / fromIntegral k
      where
        k :: Int
        k = fromIntegral i - fromIntegral j
