{-# LANGUAGE TypeApplications #-}

-- import qualified Data.Vector as V
import           Control.Monad
import           Data.Maybe
import           Numeric.Spline
import           Test.HUnit
import qualified Data.Map       as M
import qualified Data.Set       as S

main :: IO ()
main = void . runTestTT $ TestList
    [ "Sine spline" ~: splineTest
    ]

splineTest :: Assertion
splineTest = do
    expected <- map (read @Double) . lines <$> readFile "test-data/sintest.csv"
    expected @=? samples
  where
    spline :: Spline Double
    spline = fromJust
           . makeSpline
           . M.fromSet sin
           . S.fromList
           $ [0, 1, 2.5, 3.6, 5, 7, 8.1, 10]
    samples :: [Double]
    samples = sampleSpline spline . (/ 4) . fromInteger <$> [0..40]

-- test1 = (8 :: Int) ~=? (3 * 4)
--     -- TestCase $ assertEqual "4 times 2 is 8" 8 (4 * 2)
-- -- doublingBigger = TestCase $ 
-- -- halvingSmaller = TestCase "Half of 9 is 4" $ assertEqual [] 9 (9 `div` 2)

