{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE TypeApplications                         #-}
{-# LANGUAGE TypeInType                               #-}
{-# LANGUAGE TypeOperators                            #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}

import           Control.DeepSeq
import           Control.Exception
import           Criterion.Main
import           Data.Char
import           Data.Complex
import           Data.Maybe
import           GHC.TypeNats
import           Numeric.EMD
import           Text.Printf
import qualified Data.Vector       as UV
import qualified Data.Vector.Sized as V
import qualified Numeric.FFT       as FFT
import qualified System.Random.MWC as MWC

main :: IO ()
main = do
    g     <- MWC.initialize
           . UV.fromList
           . map (fromIntegral . ord)
           $ "hello world"

    test256    <- evaluate . force =<< generateData @8  g
    test1024   <- evaluate . force =<< generateData @10 g
    test4096   <- evaluate . force =<< generateData @12 g
    test16384  <- evaluate . force =<< generateData @14 g

    itest256   <- evaluate . force $ emd defaultEO test256
    itest1024  <- evaluate . force $ emd defaultEO test1024
    itest4096  <- evaluate . force $ emd defaultEO test4096
    itest16384 <- evaluate . force $ emd defaultEO test16384

    let imfs256   = length . emdIMFs $ itest256
        imfs1024  = length . emdIMFs $ itest1024
        imfs4096  = length . emdIMFs $ itest4096
        imfs16384 = length . emdIMFs $ itest16384

    defaultMainWith defaultConfig [
        bgroup "emd"
          [ bench (printf "256 (%d imfs)"   imfs256  ) $ nf (emd defaultEO) test256
          , bench (printf "1024 (%d imfs)"  imfs1024 ) $ nf (emd defaultEO) test1024
          , bench (printf "4096 (%d imfs)"  imfs4096 ) $ nf (emd defaultEO) test4096
          , bench (printf "16384 (%d imfs)" imfs16384) $ nf (emd defaultEO) test16384
          ]
      , bgroup "iemd"
          [ bench "256"   $ nf iemd itest256
          , bench "1024"  $ nf iemd itest1024
          , bench "4096"  $ nf iemd itest4096
          , bench "16384" $ nf iemd itest16384
          ]
      ]

generateData
    :: KnownNat n
    => MWC.GenIO
    -> IO (V.Vector (2^n) Double)
generateData g = fmap (fmap realPart . ifftSized) . V.generateM $ \i ->
    let i' = recip . (+ 1) . fromIntegral $ i
    in  (:+) <$> MWC.uniformR (-i', i') g
             <*> MWC.uniformR (-i', i') g
  where


ifftSized
    :: KnownNat n
    => V.Vector (2^n) (Complex Double)
    -> V.Vector (2^n) (Complex Double)
ifftSized = fromJust
          . V.fromList
          . FFT.ifft
          . V.toList
