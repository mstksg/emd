{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE TemplateHaskell                          #-}
{-# LANGUAGE TypeApplications                         #-}
{-# LANGUAGE TypeInType                               #-}
{-# LANGUAGE TypeOperators                            #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise       #-}

module Tests.EMD (
    emdTests
  ) where

import           Data.Functor.Identity
import           Data.Proxy
import           GHC.TypeNats
import           Hedgehog
import           Numeric.EMD
import           Test.Tasty
import           Tests.Util

emdTests :: TestTree
emdTests = groupTree $$(discover)

prop_iemd :: Property
prop_iemd = property $ withSize $ \(_ :: Proxy n) -> do
    xs <- forAll $ generateData @n
    tripping (CE xs) (emd @_ @_ @(2^n-1) defaultEO . getCE) (Identity . CE . iemd)
