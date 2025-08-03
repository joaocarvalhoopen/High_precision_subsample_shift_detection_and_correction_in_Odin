// Name    : High precision subsample shift detection and correction in Odin
// Date    : 2025.08.3
// Autor   : Joao Carvalho
// License : MIT Open Source

// Example output to N = 8196 sample vectors

// Begin subsample vector 1D shift detection and correction ...
//
// N : 8192
//
// Integer shift from cross-correlation peak:  -15 samples
// Calculated subsample shift :                -15.1234567889414393 samples
// True shift :                                -15.1234567890123461 samples
// Error :                                     7.0906835958339798e-11 samples
// Calculated subsample shift ( after shift ): -0.0000000001314220 samples
//
// Integer shift from cross-correlation peak:  19 samples
// Calculated subsample shift :                19.1234567888950551 samples
// True shift :                                19.1234567890123444 samples
// Error :                                     1.1728928939191974e-10 samples
// Calculated subsample shift ( after shift ): 0.0000000000209184 samples
//
// Integer shift from cross-correlation peak:  100 samples
// Calculated subsample shift :                100.1234567889132450 samples
// True shift :                                100.1234567890123515 samples
// Error :                                     9.9106500783818774e-11 samples
// Calculated subsample shift ( after shift ): 0.0000000001664375 samples
//
// Integer shift from cross-correlation peak:  -100 samples
// Calculated subsample shift :                -100.1234567890396647 samples
// True shift :                                -100.1234567890123515 samples
// Error :                                     2.7313262762618251e-11 samples
// Calculated subsample shift ( after shift ): 0.0000000000618456 samples
//
// Integer shift from cross-correlation peak:  -192 samples
// Calculated subsample shift :                8000.1234567889187019 samples
// True shift :                                8000.1234567890123799 samples
// Error :                                     9.3677954282611609e-11 samples
// Calculated subsample shift ( after shift ): 0.0000000000281943 samples
//
// Integer shift from cross-correlation peak:  192 samples
// Calculated subsample shift :                -8000.1234567890023754 samples
// True shift :                                -8000.1234567890123799 samples
// Error :                                     1.0004441719502211e-11 samples
// Calculated subsample shift ( after shift ): -0.0000000000009095 samples
//
//
// ... end subsample vector 1D shift detection and correction.
//
//
// real	0m0,122s
// user	0m0,115s
// sys	0m0,007s


package main

import "core:fmt"
import "core:math"
import "core:math/cmplx"
import "core:slice"

// A recursive implementation of the 1D Cooley-Tukey FFT algorithm.
fft_transform :: proc ( x : [ ]complex128 ) ->
                      ( res : [ ]complex128 ) {

    N := len( x )
    if N <= 1 {

        return slice.clone( x )
    }

    even_tmp := make( [ ]complex128, ( N + 1 ) / 2 )
    odd_tmp  := make( [ ]complex128, N / 2 )
    defer delete( even_tmp )
    defer delete( odd_tmp )

    for i in 0 ..= N / 2 - 1 {

        even_tmp[ i ] = x[ 2 * i ]
        odd_tmp[ i ] = x[ 2 * i  + 1 ]
    }

    if N % 2 != 0 {

        even_tmp[ N / 2 ] = x[ N - 1 ]
    }

    even := fft_transform( even_tmp )
    odd  := fft_transform( odd_tmp )
    defer delete( even )
    defer delete( odd )

    res = make( [ ]complex128, N )
    for k in 0..< N / 2 {

        angle := -2 * math.PI * f64( k ) / f64( N )
        t := cmplx.exp( complex( 0, angle ) ) * odd[ k ]
        res[ k ]         = even[ k ] + t
        res[ k + N / 2 ] = even[ k ] - t
    }
    return res
}

// A recursive implementation of the 1D Cooley-Tukey IFFT algorithm.
ifft_transform :: proc ( X : [ ]complex128 ) ->
                       ( res : [ ]complex128 ) {

    N := len( X )
    if N <= 1 {

        return slice.clone( X )
    }

    X_conj := make( [ ]complex128, N )
    for v, i in X {

        X_conj[ i ] = cmplx.conj( v )
    }
    defer delete( X_conj )

    fft_res := fft_transform( X_conj )
    defer delete( fft_res )

    res = make( [ ]complex128, N )
    for v, i in fft_res {

        res[ i ] = cmplx.conj( v ) / complex( f64( N ), 0 )
    }

    return res
}

// Calculates the cross-correlation of two 1D complex vectors using the
// Fourier transform.
//
// Cross-correlation measures the similarity between two signals as a
// function of the displacement of one relative to the other. The peak
// of the cross-correlation indicates the displacement where the signals
// are most aligned. The calculation is performed efficiently in the
// frequency domain via the formula:
//     Cross-Correlation( A, B ) = IFFT( conj( FFT( A ) ) * FFT( B ) )
//
cross_correlation :: proc ( A, B : [ ]complex128 ) ->
                          ( cross_corr : [ ]complex128 ) {

    assert( len( A ) == len( B ) )

    fft_A := fft_transform( A )
    fft_B := fft_transform( B )
    defer delete( fft_A )
    defer delete( fft_B )

    cross_power_spectrum := make( [ ]complex128, len( A ) )
    defer delete( cross_power_spectrum )
    for i in 0 ..< len( A ) {

        cross_power_spectrum[ i ] = cmplx.conj( fft_A[ i ] ) * fft_B[ i ]
    }

    cross_corr = ifft_transform( cross_power_spectrum )
    return cross_corr
}

// Performs a simple linear regression from first principles on a
// set of 2D points (x, y).
//
// This function calculates the slope (m) and y-intercept (c) of a
// line of best fit that minimizes the sum of the squared differences
// between the observed y-values and the values predicted by the line.
// It is used here to find the slope of the phase-vs-frequency plot,
// which is directly proportional to the sub-sample shift.
//
// x - A slice of f64 values for the independent variable.
// y -  A slice of f64 values for the dependent variable.
//
// return m - The slope of the regression line.
// return c - The y-intercept of the regression line.
//
linear_regression :: proc ( x, y : [ ]f64 ) ->
                          ( m, c : f64 ) {

    assert( len( x ) == len( y ) )
    N := f64( len( x ) )
    if N == 0 {

        return 0, 0
    }

    sum_x, sum_y, sum_xy, sum_x_squared : f64

    for i in 0 ..< len( x ) {

        sum_x += x[ i ]
        sum_y += y[ i ]
        sum_xy += x[ i ] * y[ i ]
        sum_x_squared += x[ i ] * x[ i ]
    }

    denominator := ( N * sum_x_squared - sum_x * sum_x )
    if denominator == 0 {

        m = 0
    } else {

        m = ( N * sum_xy - sum_x * sum_y ) / denominator
    }

    c = ( sum_y - m * sum_x ) / N

    return m, c
}

// Generates the frequency bins for a Discrete Fourier Transform ( DFT ).
//
// This function creates an array of frequencies corresponding to each
// element of an FFT output. The ordering matches standard FFT conventions:
//     the first half contains the positive frequencies ( from 0 up to the
//     Nyquist frequency ), and the second half contains the negative
//     frequencies ( from most negative back towards 0 ).
//
// n - The number of samples in the FFT.
// d - The sample spacing ( defaults to 1.0 ).
//
// return freq - A new slice of f64 values representing the frequency for each bin.
//
fft_freq :: proc ( n : int,
                   d : f64 = 1.0 ) ->
                 ( freq : [ ]f64 ) {

    val := 1.0 / ( f64( n ) * d )
    freq = make( [ ]f64, n )

    N_positive := ( n + 1 ) / 2
    N_negative := n / 2

    // Positive frequencies
    for i in 0 ..< N_positive {

        freq[ i ] = f64( i ) * val
    }

    // Negative frequencies
    for i in 0 ..< N_negative {

        freq[ N_positive + i ] = f64( -N_negative + i ) * val
    }

    return freq
}

// Corrects phase wrapping in an array of angles.
//
// The phase of a complex number is typically calculated within the
// range [ -PI, +PI ]. When the true phase moves outside this range, it
// "wraps around," creating a large jump ( e.g., from +3.1 to -3.1 ).
// This function detects these 2*PI discontinuities and adds or
// subtracts the appropriate multiple of 2*PI to make the phase
// values continuous. This is essential for performing linear
// regression on phase data.
//
// wrapped -  A slice of f64 phase values (in radians) that may be wrapped.
//
// return unwrapped - A new slice of f64 values with the phase unwrapped.
//
unwrap_phase :: proc ( wrapped : [ ]f64 ) ->
                     ( unwrapped : [ ]f64 ) {

    N := len( wrapped )
    if N == 0 {

        return nil
    }

    unwrapped = make( [ ]f64, N )
    unwrapped[ 0 ] = wrapped[ 0 ]
    offset := 0.0

    for i in 1 ..< N {

        diff := wrapped[ i ] - wrapped[ i - 1 ]
        if diff > math.PI {

            offset -= 2 * math.PI
        } else if diff < -math.PI {

            offset += 2 * math.PI
        }

        unwrapped[ i ] = wrapped[ i ] + offset
    }

    return unwrapped
}

// Calculates the sub-sample translational shift between two 1D signals.
//
// This is the core algorithm for high-precision shift detection. It
// combines cross-correlation with phase analysis to find both the
// integer and fractional parts of the shift. The steps are:
//
//   1. Compute the cross-power spectrum: `conj( FFT( A ) ) * FFT( B )`.
//
//   2. Find the integer part of the shift by locating the peak of the
//      cross-correlation ( IFFT of the cross-power spectrum ).
//
//   3. Correct the phase of the cross-power spectrum by removing the
//      linear phase component corresponding to the integer shift.
//
//   4. Unwrap the resulting phase to remove 2 * PI discontinuities.
//
//   5. Perform a linear regression on the unwrapped phase versus frequency.
//      The slope of this line is directly proportional to the sub-sample shift.
//
//   6. Combine the integer and sub-sample shifts to get the total shift.
//
// A -  The reference signal.
// B -  The shifted signal.
//
// return total_shift - The total calculated shift ( integer + fractional part ).
//
shift_detection :: proc ( A, B : [ ]complex128 ) ->
                        ( total_shift : f64 ) {

    N := len( A )
    assert( N == len( B ) )

    fft_A := fft_transform( A )
    fft_B := fft_transform( B )
    defer delete( fft_A )
    defer delete( fft_B )

    cross_power_spectrum := make( [ ]complex128, N )
    defer delete( cross_power_spectrum )
    for i in 0 ..< N {

        cross_power_spectrum[ i ] = cmplx.conj( fft_A[ i ] ) * fft_B[ i ]
    }

    cross_corr := ifft_transform( cross_power_spectrum )
    defer delete( cross_corr )

    max_abs : f64 = -1
    int_shift_idx : int
    for v, i in cross_corr {

        abs_v := cmplx.abs( v )
        if abs_v > max_abs {

            max_abs = abs_v
            int_shift_idx = i
        }
    }

    int_shift := int( int_shift_idx )
    // Handle the case where the shift is in the second half of the array
    if int_shift_idx > N / 2 {

        int_shift = int_shift_idx - N
    } else {

        int_shift = int_shift_idx
    }

    freq_vec := fft_freq( N )
    defer delete( freq_vec )

    wrapped_phase := make( [ ]f64, N )
    defer delete( wrapped_phase )

    for i in 0 ..< N {

        angle := 2 * math.PI * f64( int_shift ) * freq_vec[ i ]
        phase_correction := cmplx.exp( complex( 0, angle ) )
        corrected_val := cross_power_spectrum[ i ] * phase_correction
        wrapped_phase[ i ] = cmplx.phase( corrected_val )
    }

    // Unwrap the phase before linear regression.
    unwrapped_phase := unwrap_phase( wrapped_phase )
    defer delete( unwrapped_phase )

    m, _ := linear_regression( freq_vec, unwrapped_phase )
    subsample_shift := -m / ( 2 * math.PI )

    // The total shift must be adjusted for periodicity.
    total_shift = ( f64( int_shift ) + subsample_shift )

    // Correct for the periodic boundary
    if N > 0 {

        total_shift = math.mod( total_shift + f64( N ) / 2, f64( N ) ) - f64( N ) / 2
    }

    return total_shift
}

// Applies a precise, sub-sample shift to a signal.
//
// This function uses the Fourier Shift Theorem, which states that a
// shift in the time domain corresponds to adding a linear phase ramp
// in the frequency domain. The process is:
//
//   1. Transform the signal to the frequency domain ( FFT ).
//
//   2. Multiply each frequency component by a complex exponential
//      ( a phase shift ) whose phase is proportional to the frequency
//      and the desired shift.
//
//   3. Transform the result back to the time domain ( IFFT ).
//
// vector - The original time-domain signal.
// shift - The desired shift (can be a fractional value).
//
// return shifted_vector -  A new slice containing the shifted signal.
//
shift_vector :: proc ( vector : [ ]complex128,
                       shift  : f64 ) ->
                     ( shifted_vector : [ ]complex128 ) {

    N := len( vector )
    fft_vec := fft_transform( vector )
    defer delete( fft_vec )

    freq := fft_freq( N )
    defer delete( freq )

    shifted_fft_vec := make( [ ]complex128, N )
    defer delete( shifted_fft_vec )

    for i in 0 ..< N {

        angle := -2 * math.PI * shift * freq[ i ]
        phase_shift := cmplx.exp( complex( 0, angle ) )
        shifted_fft_vec[ i ] = fft_vec[ i ] * phase_shift
    }

    shifted_vector = ifft_transform( shifted_fft_vec )

    return shifted_vector
}


test_shift_detection_and_correction :: proc ( true_shift : f64,
                                              N          : int ) {

    A := make( [ ]complex128, N )
    defer delete( A )

    for i in 0 ..< N {
        t := f64( i )
        val := math.exp( -math.pow( t - f64( N ) / 2, 2 ) / ( 2 * math.pow( f64( N ) / 8, 2 ) ) )
        A[ i ] = complex( val, 0 )
    }

    B := shift_vector( A, true_shift )
    defer delete( B )

    cross_corr := cross_correlation( A, B )
    defer delete( cross_corr )

    max_abs : f64 = -1
    int_shift_idx : int
    for v, i in cross_corr {

        abs_v := cmplx.abs( v )
        if abs_v > max_abs {

            max_abs = abs_v
            int_shift_idx = i
        }
    }

    integer_shift := int( int_shift_idx )
    if integer_shift > N / 2 {

        integer_shift -= N
    }

    fmt.printfln( "Integer shift from cross-correlation peak:  %d samples", integer_shift )

    total_shift := shift_detection( A, B )

    // We need to adjust the final shift to be in the correct periodic range for comparison
    final_shift_for_comparison := total_shift

    if true_shift > f64( N ) / 2 && total_shift < 0 {

        final_shift_for_comparison += f64( N )
    } else if true_shift < -f64( N ) / 2 && total_shift > 0 {

        final_shift_for_comparison -= f64( N )
    }

    fmt.printfln( "Calculated subsample shift :                %.16f samples", final_shift_for_comparison )
    fmt.printfln( "True shift :                                %.16f samples", true_shift )
    fmt.printfln( "Error :                                     %.16e samples", math.abs( final_shift_for_comparison - true_shift ) )

    C := shift_vector( B, -total_shift )
    defer delete( C )

    total_shift_last := shift_detection( A, C )
    fmt.printfln( "Calculated subsample shift ( after shift ): %.16f samples\n", total_shift_last )
}


main :: proc ( ) {

    fmt.println( "\nBegin subsample vector 1D shift detection and correction ...\n" )

    N :: 8192 // 4096

    fmt.printfln( "N : %d\n", N )

    true_shift : f64 = -15.12345678901234567890
    test_shift_detection_and_correction( true_shift, N )

    true_shift = +19.12345678901234567890
    test_shift_detection_and_correction( true_shift, N )

    true_shift = +100.12345678901234567890
    test_shift_detection_and_correction( true_shift, N )

    true_shift = -100.12345678901234567890
    test_shift_detection_and_correction( true_shift, N )

    true_shift = +4000.12345678901234567890
    test_shift_detection_and_correction( true_shift, N )

    true_shift = -4000.12345678901234567890
    test_shift_detection_and_correction( true_shift, N )

    fmt.println( "\n ... end subsample vector 1D shift detection and correction.\n" )
}
