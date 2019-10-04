Changelog
=========

Version 0.1.8.0
---------------

*October 4, 2019*

<https://github.com/mstksg/emd/releases/tag/v0.1.8.0>

*   Add `meanMarginal`
*   Fix `degreeOfStationarity` for divide-by-zero errors.
*   Add `foldFreq` for generalized folding on `HHT`, and rewrote other
    functions in terms of it.
*   Drop support for GHC 8.2 and lower.

*   *0.1.8.1*: Exported `marginal` again; it was unexported by mistake.

Version 0.1.7.0
---------------

*September 24, 2019*

<https://github.com/mstksg/emd/releases/tag/v0.1.7.0>

*   Rewrite `hilbert` using the *fft* library, matching the matlab
    implementation.  This means that the library now depends on *fftw*.

Version 0.1.6.0
---------------

*September 24, 2019*

<https://github.com/mstksg/emd/releases/tag/v0.1.6.0>

*   Add `hilbertPhase` to *Numeric.HHT*.

Version 0.1.5.1
---------------

*September 3, 2019*

<https://github.com/mstksg/emd/releases/tag/v0.1.5.1>

*   Remove dependency on *pure-fft*, using *statistics* instead.

Version 0.1.5.0
---------------

*August 31, 2019*

<https://github.com/mstksg/emd/releases/tag/v0.1.5.0>

*   Add `NFData` instance for `EMD`, `HHTLine`, and `HTT`
*   Add `iemd`, inverting `emd`.

Version 0.1.4.0
---------------

*August 20, 2018*

<https://github.com/mstksg/emd/releases/tag/v0.1.4.0>

*   `hhtSparseSpectrum` added to *Numeric.HHT* module, for an alternate sparser
    representation of the Hilbert Spectrum.
*   `hhtDenseSpectrum` also added to *Numeric.HHT*, for an alternative denser
    representation.
*   `expectedFrequency` added to *Numeric.HHT* module, to calculate weighted
    average of frequency contributions at each step in time.
*   `dominantFrequency` also added to *Numeric.HHT* to calculate strongest
    frequency at each step in time.

Version 0.1.3.0
---------------

*August 15, 2018*

<https://github.com/mstksg/emd/releases/tag/v0.1.3.0>

*   `Default` instance for `SiftCondition` and `EMDOpts`, as a useful
    alternative to `defaultEO` and `defaultSC` for those who prefer it.
*   `Binary` instances for `EMD`, `HHT`, and related data types.  These are
    based on `Binary` instance for `v a`, so the user must bring the orphan
    instance of their choice into scope.  Not sure if this is the best way to
    do this.

Version 0.1.2.1
---------------

*July 27, 2018*

<https://github.com/mstksg/emd/releases/tag/v0.1.2.1>

*   *BUG FIX* Fixed behavior of frequency wrapping to wrap between 0 and 1,
    instead of 0.5, as claimed!

Version 0.1.2.0
---------------

*July 27, 2018*

<https://github.com/mstksg/emd/releases/tag/v0.1.2.0>

*   Actually implemented the Hilbert-Huang transform
*   Allowed for other border handling behaviors during EMD
*   Changed default stopping conditions for sifting process
*   Added clamped spline end conditions.
*   Removed unsized interface
*   Sifting will now throw runtime exception for singular splining matrices,
    instead of treating the result as a residual.  This might change in the
    future.

Version 0.1.1.0
---------------

*July 25, 2018*

<https://github.com/mstksg/emd/releases/tag/v0.1.1.0>

*   Unsized interface added.

Version 0.1.0.0
---------------

*July 25, 2018*

<https://github.com/mstksg/emd/releases/tag/v0.1.0.0>

*   Initial release
