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

prop_iemd_default :: Property
prop_iemd_default = iemdProp defaultEO

iemdProp :: EMDOpts Double -> Property
iemdProp eo = property $ withSize $ \(_ :: Proxy n) -> do
    xs <- forAll $ generateData @n
    tripping (CE xs) (emd @_ @_ @(2^n-1) eo . getCE) (Identity . CE . iemd)
