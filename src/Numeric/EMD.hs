{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE LambdaCase   #-}

module Numeric.EMD (
    locals
  ) where

import           Control.Monad.Trans.State
import           Control.Monad.Trans.Writer
import           Data.Finite
import           Data.Monoid
import           Data.Ord
import           GHC.TypeNats
import qualified Data.Map                   as M
import qualified Data.Vector.Generic        as VG
import qualified Data.Vector.Generic.Sized  as SVG

locals
    :: (VG.Vector v a, KnownNat n, Num a, Ord a)
    => SVG.Vector v n a
    -> (M.Map (Finite n) a, M.Map (Finite n) a)
locals = snd . SVG.ifoldl' (uncurry go) (Nothing, (M.empty, M.empty))
  where
    go = \case
      Nothing            -> \mms i x ->  (Just (x, EQ, i), mms)
      Just (!y, !d, !i') -> \mms@(!mins, !maxs) !i !x ->
        let d'  = (x - y) `compare` 0
            mms' = case (d, d') of
              (LT, LT) -> mms
              (LT, EQ) -> (M.insert i' y mins, maxs)
              (LT, GT) -> (M.insert i' y mins, maxs)
              (EQ, _ ) -> mms
              (GT, LT) -> (mins, M.insert i' y maxs)
              (GT, EQ) -> (mins, M.insert i' y maxs)
              (GT, GT) -> mms
        in  (Just (x, d', i), mms')


