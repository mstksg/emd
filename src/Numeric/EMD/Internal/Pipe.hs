{-# LANGUAGE DeriveFunctor       #-}
{-# LANGUAGE LambdaCase          #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Numeric.EMD.Internal.Pipe (
    Pipe
  , (.|)
  , runPipe
  , awaitSurely
  , repeatM
  , dropP
  , ZipSink(..)
  -- , sinkList
  -- , interleaveP
  -- , yield
  -- , await
  -- , awaitForever
  -- , sourceList
  -- , unfoldP, unfoldPForever, sourceList, iterateP
  -- , mapP
  -- , mapMP
  -- , takeP
  -- , sinkList
  ) where

import           Control.Applicative
import           Control.Monad
import           Data.Void
import           Control.Monad.Trans.Class

-- | Pipe
--
-- *  @i@: Input
-- *  @o@: Output
-- *  @u@: Upstream result
-- *  @m@: Monad
-- *  @r@: Result
--
-- Some specializations:
--
-- *  If @i@ is '()', we have a producer.  If @r@ is 'Void', it will always
--    produce.
-- *  If @o@ is 'Void', we have a consumer.
-- *  If @u@ is 'Void', the pipe will never stop producing.
--
-- TODO: CPS this maybe
data Pipe i o u m a =
      PAwait (i -> Pipe i o u m a) (u -> Pipe i o u m a)
    | PYield o (Pipe i o u m a)
    | PAct (m (Pipe i o u m a))
    | PDone a
  deriving Functor

(.|) :: Functor m => Pipe a b u m v -> Pipe b c v m r -> Pipe a c u m r
p .| q = case q of
    PAwait f g -> case p of
      PAwait f' g' -> PAwait (\x -> f' x .| q) (\x -> g' x .| q)
      PYield x' y' -> y' .| f x'
      PAct   x'    -> PAct $ (.| q) <$> x'
      PDone  x'    -> p .| g x'
    PYield x y -> PYield x (p .| y)
    PAct   x   -> PAct $ (p .|) <$> x
    PDone  x   -> PDone x
infixr 2 .|

instance Functor m => Applicative (Pipe i o u m) where
    pure = PDone
    (<*>) = \case
      PAwait f g -> \q -> PAwait ((<*> q) . f) ((<*> q) . g)
      PYield x y -> \q -> PYield x (y <*> q)
      PAct   x   -> \q -> PAct ((<*> q) <$> x)
      PDone  f   -> fmap f
    (*>) = \case
      PAwait f g -> \q -> PAwait ((*> q) . f) ((*> q) . g)
      PYield x y -> \q -> PYield x (y *> q)
      PAct   x   -> \q -> PAct ((*> q) <$> x)
      PDone  _   -> id

instance Functor m => Monad (Pipe i o u m) where
    return = PDone
    (>>=) = \case
      PAwait f g -> \q -> PAwait (f >=> q) (g >=> q)
      PYield x y -> \q -> PYield x (y >>= q)
      PAct   x   -> \q -> PAct ((>>= q) <$> x)
      PDone  x   -> ($ x)
    (>>)  = (*>)

instance MonadTrans (Pipe i o u) where
    lift = liftP

liftP :: Functor m => m a -> Pipe i o u m a
liftP = PAct . fmap PDone

runPipe :: Monad m => Pipe () Void u m a -> m a
runPipe = \case
    PAwait f _ -> runPipe $ f ()
    PYield x _ -> absurd x
    PAct   x   -> runPipe =<< x
    PDone  x   -> pure x

yield :: o -> Pipe i o u m ()
yield x = PYield x (PDone ())

awaitEither :: Pipe i o u m (Either i u)
awaitEither = PAwait (PDone . Left) (PDone . Right)

await :: Pipe i o u m (Maybe i)
await = PAwait (PDone . Just) (\_ -> PDone Nothing)

awaitSurely :: Pipe i o Void m i
awaitSurely = PAwait PDone absurd

-- unfoldP :: (b -> Maybe (a, b)) -> b -> Pipe i a u ()
-- unfoldP f = go
--   where
--     go z = case f z of
--       Nothing      -> PDone ()
--       Just (x, z') -> PYield x (go z')

-- unfoldPForever :: (b -> (a, b)) -> b -> Pipe i a u r
-- unfoldPForever f = go
--   where
--     go z = PYield x (go z')
--       where
--         (x, z') = f z

-- iterateP :: (a -> a) -> a -> Pipe i a u r
-- iterateP f = unfoldPForever (join (,) . f)

sourceList :: Foldable t => t a -> Pipe i a u m ()
sourceList = foldr PYield (PDone ())

repeatM :: Functor m => m o -> Pipe i o u m u
repeatM x = go
  where
    go = do
      yield =<< liftP x
      go

awaitForever :: Functor m => (i -> Pipe i o u m a) -> Pipe i o u m u
awaitForever f = go
  where
    go = PAwait (\x -> f x *> go) PDone

mapP :: Functor m => (a -> b) -> Pipe a b u m u
mapP f = awaitForever (yield . f)

mapMP :: Monad m => (a -> m b) -> Pipe a b u m u
mapMP f = awaitForever ((yield =<<) . lift . f)

dropP :: Functor m => Int -> Pipe i o u m ()
dropP n = replicateM_ n await
    -- awaitForever yield

takeP :: Functor m => Int -> Pipe i i u m ()
takeP n = replicateM_ n $ mapM_ yield =<< await

foldrP :: Functor m =>(a -> b -> b) -> b -> Pipe a Void u m b
foldrP f z = go
  where
    go = await >>= \case
      Nothing -> pure z
      Just x  -> f x <$> go

sinkList :: Functor m => Pipe i Void u m [i]
sinkList = foldrP (:) []

-- interleaveP :: Pipe i o u () -> Pipe i o u () -> Pipe i o u ()
-- interleaveP p = case p of
--     PAwait f g -> \case
--       PAwait f' g' -> PAwait (interleaveP <$> f <*> f') (interleaveP <$> g <*> g')
--       PYield x' y' -> PYield x' $ interleaveP p y'
--       PDone _      -> p
--     PYield x y -> \case
--       q@(PAwait _ _) -> PYield x $ interleaveP y q
--       PYield x' y'   -> PYield x . PYield x' $ interleaveP y y'
--       PDone _        -> p
--     PDone   _  -> id

newtype ZipSink i u m a = ZipSink { getZipSink :: Pipe i Void u m a }
  deriving Functor

zipSink
    :: Functor m
    => Pipe i Void u m (a -> b)
    -> Pipe i Void u m a
    -> Pipe i Void u m b
zipSink p q = case p of
    PAwait f g -> case q of
      PAwait f' g' -> PAwait (zipSink <$> f <*> f') (zipSink <$> g <*> g')
      PYield x' _  -> absurd x'
      PAct   x'    -> PAct (zipSink p <$> x')
      PDone  _     -> PAwait ((`zipSink` q) . f) ((`zipSink` q) . g)
    PYield x _ -> absurd x
    PAct   x   -> PAct ((`zipSink` q) <$> x)
    PDone  x   -> case q of
      PAwait f' g' -> PAwait (zipSink p . f') (zipSink p . g')
      PYield x' _  -> absurd x'
      PAct   x'    -> PAct (zipSink p <$> x')
      PDone  x'    -> PDone (x x')

altSink
    :: Functor m
    => Pipe i Void u m a
    -> Pipe i Void u m a
    -> Pipe i Void u m a
altSink p q = case p of
    PAwait f g -> case q of
      PAwait f' g' -> PAwait (altSink <$> f <*> f') (altSink <$> g <*> g')
      PYield x' _  -> absurd x'
      PAct   x'    -> PAct (altSink p <$> x')
      PDone  x'    -> PDone x'
    PYield x _ -> absurd x
    PAct   x   -> case q of
      PDone  x'   -> PDone x'
      _           -> PAct ((`altSink` q) <$> x)
    PDone  x   -> PDone x

-- | '<*>' = distribute input to all, and return result when they finish
--
-- 'pure' = immediately finish
instance Functor m => Applicative (ZipSink i u m) where
    pure = ZipSink . pure
    ZipSink p <*> ZipSink q = ZipSink $ zipSink p q

-- | '<|>' = distribute input to all, and return the first result that
-- finishes
--
-- 'empty' = never finish
instance Functor m => Alternative (ZipSink i u m) where
    empty = ZipSink go
      where
        go = PAwait (const go) (const go)
    ZipSink p <|> ZipSink q = ZipSink $ altSink p q
