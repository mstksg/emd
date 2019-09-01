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

    defaultMainWith defaultConfig [
        bgroup "emd"
          [ bench "256"   $ nf (emd defaultEO) test256
          , bench "1024"  $ nf (emd defaultEO) test1024
          , bench "4096"  $ nf (emd defaultEO) test4096
          , bench "16384" $ nf (emd defaultEO) test16384
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
    let i' = sqrt . recip . fromIntegral $ i
    in  (:+) <$> MWC.uniformR (-i', i') g
             <*> MWC.uniformR (-i', i') g

ifftSized
    :: KnownNat n
    => V.Vector (2^n) (Complex Double)
    -> V.Vector (2^n) (Complex Double)
ifftSized = fromJust
          . V.fromList
          . FFT.ifft
          . V.toList
