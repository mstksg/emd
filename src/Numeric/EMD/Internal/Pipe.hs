{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE LambdaCase    #-}

module Numeric.EMD.Internal.Pipe where

-- module Numeric.EMD.Internal.Pipe (
--   ) where

import           Control.Monad
import           Data.Void

data Pipe i o a =
      PAwait (Maybe i -> Pipe i o a)
    | PYield o (Pipe i o a)
    | PDone a
  deriving Functor

instance Applicative (Pipe i o) where
    pure = PDone
    (<*>) = \case
      PAwait f   -> \q -> PAwait   $ \x -> f x <*> q
      PYield x y -> \q -> PYield x (y <*> q)
      PDone  f   -> fmap f
    (*>) = \case
      PAwait f   -> \q -> PAwait   $ \x -> f x *> q
      PYield x y -> \q -> PYield x (y *> q)
      PDone  _   -> id

instance Monad (Pipe i o) where
    return = PDone
    (>>=) = \case
      PAwait f   -> \q -> PAwait   $ f >=> q
      PYield x y -> \q -> PYield x (y >>= q)
      PDone  x -> ($ x)
    (>>)  = (*>)

runPipe :: Pipe () Void a -> a
runPipe = \case
    PAwait f   -> runPipe $ f (Just ())
    PYield x _ -> absurd x
    PDone  x   -> x

yield :: o -> Pipe i o ()
yield x = PYield x (pure ())

await :: Pipe i o (Maybe i)
await = PAwait PDone

sourceList :: Foldable t => t a -> Pipe i a ()
sourceList = foldr PYield (PDone ())

awaitForever :: (i -> Pipe i o a) -> Pipe i o ()
awaitForever f = go
  where
    go = PAwait $ \case
      Nothing -> PDone ()
      Just x  -> f x *> go

mapP :: (a -> b) -> Pipe a b ()
mapP f = awaitForever (yield . f)

foldrP :: (a -> b -> b) -> b -> Pipe a Void b
foldrP f z = go
  where
    go = await >>= \case
      Nothing -> pure z
      Just x  -> f x <$> go

sinkList :: Pipe i Void [i]
sinkList = foldrP (:) []

headP :: Pipe i Void (Maybe i)
headP = await

(.|) :: Pipe a b () -> Pipe b c r -> Pipe a c r
(.|) = \case
    PAwait f   -> \q -> PAwait $ \x -> f x .| q
    PYield x y -> \case
      PAwait f     -> y .| f (Just x)
      PYield x' y' -> PYield x' (y .| y')
      PDone  r     -> PDone r
    PDone  _   -> \case
      PAwait f   -> PAwait $ \_ -> PDone () .| f Nothing
      PYield x y -> PYield x (PDone () .| y)
      PDone  r   -> PDone r
infixr 2 .|


