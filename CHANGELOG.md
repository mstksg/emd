Changelog
=========

Version 0.2.0.0
---------------

*August 13, 2018*

<https://github.com/mstksg/emd/releases/tag/v0.2.0.0>

*   Rewrite API to be "more type-safe".  `EMD` and `HHT` now contain the number
    of IMFs in their type, and `emd` (and `emdTrace` and `emd'`) and `hht` now
    return existentially quantified types in a continuation.

    The only real benefit is that `hhtEmd` now guarantees that it preserves the
    number of IMFs.

    You can use `someEmd`, `someHht` with the `SomeEMD` and `SomeHHT` wrapper
    types to recover the original API.


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
