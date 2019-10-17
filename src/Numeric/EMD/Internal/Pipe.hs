{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE LambdaCase    #-}

module Numeric.EMD.Internal.Pipe (

    Pipe(..)
  , (.|)
  , runPipe
  , yield, await, awaitSurely
  , unfoldP, unfoldPForever, sourceList, iterateP
  , mapP
  , sinkList
  , interleaveP
  ) where

import           Control.Monad
import           Data.Void

-- | Pipe
--
-- *  @i@: Input
-- *  @o@: Output
-- *  @u@: Upstream result
-- *  @r@: Result
--
-- Some specializations:
--
-- *  If @i@ is '()', we have a producer.  If @r@ is 'Void', it will always
--    produce.
-- *  If @o@ is 'Void', we have a consumer.
-- *  If @u@ is 'Void', the pipe will never stop producing.
--
data Pipe i o u a =
      PAwait (i -> Pipe i o u a) (u -> Pipe i o u a)
    | PYield o (Pipe i o u a)
    | PDone a
  deriving Functor

(.|) :: Pipe a b u v -> Pipe b c v r -> Pipe a c u r
(.|) = \case
    PAwait f g -> \q -> PAwait (\x -> f x .| q) (\x -> g x .| q)
    PYield x y -> \case
      PAwait f _   -> y .| f x
      PYield x' y' -> PYield x' (y .| y')
      PDone  r     -> PDone r
    r@(PDone  r') -> \case
      PAwait _ g -> r .| g r'
      PYield x y -> PYield x (r .| y)
      PDone  s   -> PDone s
infixr 2 .|

instance Applicative (Pipe i o u) where
    pure = PDone
    (<*>) = \case
      PAwait f g -> \q -> PAwait (\x -> f x <*> q) (\x -> g x <*> q)
      PYield x y -> \q -> PYield x (y <*> q)
      PDone  f   -> fmap f
    (*>) = \case
      PAwait f g -> \q -> PAwait (\x -> f x *> q) (\x -> g x *> q)
      PYield x y -> \q -> PYield x (y *> q)
      PDone  _   -> id

instance Monad (Pipe i o u) where
    return = PDone
    (>>=) = \case
      PAwait f g -> \q -> PAwait (f >=> q) (g >=> q)
      PYield x y -> \q -> PYield x (y >>= q)
      PDone  x -> ($ x)
    (>>)  = (*>)

runPipe :: Pipe () Void u a -> a
runPipe = \case
    PAwait f _ -> runPipe $ f ()
    PYield x _ -> absurd x
    PDone  x   -> x

yield :: o -> Pipe i o u ()
yield x = PYield x (pure ())

await :: Pipe i o u (Maybe i)
await = PAwait (PDone . Just) (\_ -> PDone Nothing)

awaitSurely :: Pipe i o Void i
awaitSurely = PAwait PDone absurd

unfoldP :: (b -> Maybe (a, b)) -> b -> Pipe i a u ()
unfoldP f = go
  where
    go z = case f z of
      Nothing      -> PDone ()
      Just (x, z') -> PYield x (go z')

unfoldPForever :: (b -> (a, b)) -> b -> Pipe i a u r
unfoldPForever f = go
  where
    go z = PYield x (go z')
      where
        (x, z') = f z

iterateP :: (a -> a) -> a -> Pipe i a u r
iterateP f = unfoldPForever (join (,) . f)

sourceList :: Foldable t => t a -> Pipe i a u ()
sourceList = foldr PYield (PDone ())

awaitForever :: (i -> Pipe i o u a) -> Pipe i o u ()
awaitForever f = go
  where
    go = PAwait (\x -> f x *> go) (\_ -> PDone ())

mapP :: (a -> b) -> Pipe a b u ()
mapP f = awaitForever (yield . f)

foldrP :: (a -> b -> b) -> b -> Pipe a Void u b
foldrP f z = go
  where
    go = await >>= \case
      Nothing -> pure z
      Just x  -> f x <$> go

sinkList :: Pipe i Void u [i]
sinkList = foldrP (:) []

interleaveP :: Pipe i o u () -> Pipe i o u () -> Pipe i o u ()
interleaveP p = case p of
    PAwait f g -> \case
      PAwait f' g' -> PAwait (interleaveP <$> f <*> f') (interleaveP <$> g <*> g')
      PYield x' y' -> PYield x' $ interleaveP p y'
      PDone _      -> p
    PYield x y -> \case
      q@(PAwait _ _) -> PYield x $ interleaveP y q
      PYield x' y'   -> PYield x . PYield x' $ interleaveP y y'
      PDone _        -> p
    PDone   _  -> id
