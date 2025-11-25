/*
 * Copyright ConsenSys AG.
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
 * an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 *
 * SPDX-License-Identifier: Apache-2.0
 */
package org.hyperledger.besu.evm;

import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;

import com.google.common.annotations.VisibleForTesting;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * 256-bits wide unsigned integer class.
 *
 * <p>This class is an optimised version of BigInteger for fixed width 256-bits integers.
 */
public final class UInt256 {
  // region Internals
  // --------------------------------------------------------------------------
  // UInt256 is a big-endian up to 256-bits integer.
  // Internally, it is represented with fixed-size int/long limbs in little-endian order.
  // Length is used to optimise algorithms, skipping leading zeroes.
  // Nonetheless, 256bits are always allocated and initialised to zeroes.

  /** Fixed size in bytes. */
  public static final int BYTESIZE = 32;

  /** Fixed size in bits. */
  public static final int BITSIZE = 256;

  // Fixed number of limbs or digits
  private static final int N_LIMBS = 8;
  // Fixed number of bytes per limb.
  private static final int N_BYTES_PER_LIMB = 4;
  // Fixed number of bits per limb.
  private static final int N_BITS_PER_LIMB = 32;
  // Mask for long values
  private static final long MASK_L = 0xFFFFFFFFL;

  // SIMD Vector Species
  private static VectorSpecies<Integer> SPECIES = IntVector.SPECIES_256;
  private static VectorSpecies<Long> SPECIES_L = LongVector.SPECIES_256;

  // Vectors utils
  private static LongVector ALLSET_VEC = LongVector.broadcast(SPECIES_L, -1L); // 0xFFFFFFFFFFFFFFFF

  // Internal state
  private final int[] limbs;
  private final int length;

  @VisibleForTesting
  int[] limbs() {
    return limbs;
  }

  // --------------------------------------------------------------------------
  // endregion

  /** The constant 0. */
  public static final UInt256 ZERO = new UInt256(new int[] {0, 0, 0, 0, 0, 0, 0, 0}, 0);

  // region Constructors
  // --------------------------------------------------------------------------

  UInt256(final int[] limbs, final int length) {
    // Unchecked length: assumes limbs have length == N_LIMBS
    this.limbs = limbs;
    this.length = length;
  }

  UInt256(final int[] limbs) {
    this(limbs, N_LIMBS);
  }

  /**
   * Instantiates a new UInt256 from byte array.
   *
   * @param bytes raw bytes in BigEndian order.
   * @return Big-endian UInt256 represented by the bytes.
   */
  public static UInt256 fromBytesBE(final byte[] bytes) {
    int nLimbs = (bytes.length + N_BYTES_PER_LIMB - 1) / N_BYTES_PER_LIMB;
    int nBytes = nLimbs * N_BYTES_PER_LIMB;
    byte[] padded = new byte[nBytes];
    System.arraycopy(bytes, 0, padded, nBytes - bytes.length, bytes.length);
    ByteBuffer buf = ByteBuffer.wrap(padded).order(ByteOrder.BIG_ENDIAN);
    int[] limbs = new int[N_LIMBS];
    for (int i = nLimbs - 1; i >= 0; i--) {
      limbs[i] = buf.getInt(); // reverse order for little-endian limbs
    }
    return new UInt256(limbs, nLimbs);
  }

  /**
   * Instantiates a new UInt256 from an int.
   *
   * @param value int value to convert to UInt256.
   * @return The UInt256 equivalent of value.
   */
  public static UInt256 fromInt(final int value) {
    if (value == 0) return ZERO;
    int[] limbs = new int[N_LIMBS];
    limbs[0] = value;
    return new UInt256(limbs, 1);
  }

  /**
   * Instantiates a new UInt256 from a long.
   *
   * @param value long value to convert to UInt256.
   * @return The UInt256 equivalent of value.
   */
  public static UInt256 fromLong(final long value) {
    if (value == 0) return ZERO;
    int[] limbs = new int[N_LIMBS];
    limbs[0] = (int) value;
    limbs[1] = (int) (value >>> 32);
    return new UInt256(limbs, 2);
  }

  /**
   * Instantiates a new UInt256 from an int array.
   *
   * <p>The array is interpreted in little-endian order. It is either padded with 0s or truncated if
   * necessary.
   *
   * @param arr int array of limbs.
   * @return The UInt256 equivalent of value.
   */
  public static UInt256 fromArray(final int[] arr) {
    int[] limbs = new int[N_LIMBS];
    int len = Math.min(N_LIMBS, arr.length);
    System.arraycopy(arr, 0, limbs, 0, len);
    return new UInt256(limbs, len);
  }

  // --------------------------------------------------------------------------
  // endregion

  // region Conversions
  // --------------------------------------------------------------------------
  /**
   * Convert to int.
   *
   * @return Value truncated to an int, possibly lossy.
   */
  public int intValue() {
    return limbs[0];
  }

  /**
   * Convert to long.
   *
   * @return Value truncated to a long, possibly lossy.
   */
  public long longValue() {
    return (limbs[0] & MASK_L) | ((limbs[1] & MASK_L) << 32);
  }

  /**
   * Convert to BigEndian byte array.
   *
   * @return Big-endian ordered bytes for this UInt256 value.
   */
  public byte[] toBytesBE() {
    ByteBuffer buf = ByteBuffer.allocate(BYTESIZE).order(ByteOrder.BIG_ENDIAN);
    for (int i = N_LIMBS - 1; i >= 0; i--) {
      buf.putInt(limbs[i]);
    }
    return buf.array();
  }

  /**
   * Convert to BigInteger.
   *
   * @return BigInteger representing the integer.
   */
  public BigInteger toBigInteger() {
    return new BigInteger(1, toBytesBE());
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder("0x");
    for (byte b : toBytesBE()) {
      sb.append(String.format("%02x", b));
    }
    return sb.toString();
  }

  private long[] toLongLimbs() {
    int n_longs = (N_LIMBS + 1) / 2;
    long[] result = new long[n_longs];
    for (int i = 0; i < n_longs; i++) {
      result[i] = (limbs[i * 2] & MASK_L) | ((limbs[i * 2 + 1] & MASK_L) << 32);
    }
    return result;
  }

  // --------------------------------------------------------------------------
  // endregion

  // region Comparisons
  // --------------------------------------------------------------------------

  /**
   * Is the value 0 ?
   *
   * @return true if this UInt256 value is 0.
   */
  public boolean isZero() {
    return (limbs[0] | limbs[1] | limbs[2] | limbs[3] | limbs[4] | limbs[5] | limbs[6] | limbs[7])
        == 0;
  }

  /**
   * Is the value 0 ? (SIMD-accelerated)
   *
   * @return true if this UInt256 value is 0.
   */
  public boolean isZeroVec() {
    IntVector v = IntVector.fromArray(SPECIES, limbs, 0);
    return v.reduceLanes(VectorOperators.OR) == 0;
  }

  /**
   * Is the value 0 using intrinsic vectorized comparison.
   *
   * <p>This method uses Arrays.mismatch to efficiently check if
   * all limbs are zero by comparing against a static zero array.
   *
   * @return true if this UInt256 value is 0.
   */
  public boolean isZeroIntrinsic() {
    int mismatch = Arrays.mismatch(limbs, ZERO.limbs);
    return mismatch == -1;
  }

  /**
   * Compares two UInt256.
   *
   * @param a left UInt256
   * @param b right UInt256
   * @return 0 if a == b, negative if a &lt; b and positive if a &gt; b.
   */
  public static int compare(final UInt256 a, final UInt256 b) {
    int comp;
    for (int i = N_LIMBS - 1; i >= 0; i--) {
      comp = Integer.compareUnsigned(a.limbs[i], b.limbs[i]);
      if (comp != 0) return comp;
    }
    return 0;
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) return true;
    if (!(obj instanceof UInt256)) return false;
    UInt256 other = (UInt256) obj;

    int xor =
        (this.limbs[0] ^ other.limbs[0])
            | (this.limbs[1] ^ other.limbs[1])
            | (this.limbs[2] ^ other.limbs[2])
            | (this.limbs[3] ^ other.limbs[3])
            | (this.limbs[4] ^ other.limbs[4])
            | (this.limbs[5] ^ other.limbs[5])
            | (this.limbs[6] ^ other.limbs[6])
            | (this.limbs[7] ^ other.limbs[7]);
    return xor == 0;
  }

  public boolean equalsVec(final Object obj) {
    if (this == obj) return true;
    if (!(obj instanceof UInt256)) return false;
    UInt256 other = (UInt256) obj;

    // SIMD-accelerated: XOR all limbs vectorized, then reduce with OR
    IntVector thisVec = IntVector.fromArray(SPECIES, this.limbs, 0);
    IntVector otherVec = IntVector.fromArray(SPECIES, other.limbs, 0);
    IntVector xorVec = thisVec.lanewise(VectorOperators.XOR, otherVec);
    return xorVec.reduceLanes(VectorOperators.OR) == 0;
  }

  @Override
  public int hashCode() {
    int h = 1;
    for (int i = 0; i < N_LIMBS; i++) {
      h = 31 * h + limbs[i];
    }
    return h;
  }

  // --------------------------------------------------------------------------
  // endregion

  // region Arithmetic Operations
  // --------------------------------------------------------------------------

  /**
   * Unsigned modulo reduction.
   *
   * @param modulus The modulus of the reduction
   * @return The remainder modulo {@code modulus}.
   */
  public UInt256 mod(final UInt256 modulus) {
    if (this.isZero() || modulus.isZero()) return ZERO;
    return new UInt256(knuthRemainder(this.limbs, modulus.limbs), modulus.length);
  }

  /**
   * Signed modulo reduction.
   *
   * <p>In signed modulo reduction, integers are interpretated as fixed 256 bits width two's
   * complement signed integers.
   *
   * @param modulus The modulus of the reduction
   * @return The remainder modulo {@code modulus}.
   */
  public UInt256 signedMod(final UInt256 modulus) {
    if (this.isZero() || modulus.isZero()) return ZERO;
    int[] x = new int[N_LIMBS];
    int[] y = new int[N_LIMBS];
    absInto(x, this.limbs, N_LIMBS);
    absInto(y, modulus.limbs, N_LIMBS);
    int[] r = knuthRemainder(x, y);
    if (isNeg(this.limbs, N_LIMBS)) {
      negate(r, N_LIMBS);
      return new UInt256(r);
    }
    return new UInt256(r, modulus.length);
  }

  /**
   * Modular addition.
   *
   * @param other The integer to add to this.
   * @param modulus The modulus of the reduction.
   * @return This integer this + other (mod modulus).
   */
  public UInt256 addMod(final UInt256 other, final UInt256 modulus) {
    if (modulus.isZero()) return ZERO;
    int[] sum = addWithCarry(this.limbs, this.length, other.limbs, other.length);
    int[] rem = knuthRemainder(sum, modulus.limbs);
    return new UInt256(rem, modulus.length);
  }

  /**
   * Modular multiplication.
   *
   * @param other The integer to add to this.
   * @param modulus The modulus of the reduction.
   * @return This integer this + other (mod modulus).
   */
  public UInt256 mulMod(final UInt256 other, final UInt256 modulus) {
    if (this.isZero() || other.isZero() || modulus.isZero()) return ZERO;
    int[] result = addMul(this.limbs, this.length, other.limbs, other.length);
    result = knuthRemainder(result, modulus.limbs);
    return new UInt256(result, modulus.length);
  }

  // --------------------------------------------------------------------------
  // endregion

  // region Support (private) Algorithms
  // --------------------------------------------------------------------------
  private static int nSetLimbs(final int[] x) {
    int offset = x.length - 1;
    while ((offset >= 0) && (x[offset] == 0)) offset--;
    return offset + 1;
  }

  private static int compareLimbs(final int[] a, final int aLen, final int[] b, final int bLen) {
    int cmp;
    if (aLen > bLen) {
      for (int i = aLen - 1; i >= bLen; i--) {
        cmp = Integer.compareUnsigned(a[i], 0);
        if (cmp != 0) return cmp;
      }
    } else if (aLen < bLen) {
      for (int i = bLen - 1; i >= aLen; i--) {
        cmp = Integer.compareUnsigned(0, b[i]);
        if (cmp != 0) return cmp;
      }
    }
    for (int i = Math.min(aLen, bLen) - 1; i >= 0; i--) {
      cmp = Integer.compareUnsigned(a[i], b[i]);
      if (cmp != 0) return cmp;
    }
    return 0;
  }

  private static boolean isNeg(final int[] x, final int xLen) {
    return x[xLen - 1] < 0;
  }

  private static void negate(final int[] x, final int xLen) {
    int carry = 1;
    for (int i = 0; i < xLen; i++) {
      x[i] = ~x[i] + carry;
      carry = (x[i] == 0 && carry == 1) ? 1 : 0;
    }
  }

  private static void negateVec(final int[] x, final int xLen) {
    if (xLen == N_LIMBS) {
      // SIMD-accelerated path for full-width negation
      // Step 1: Bitwise NOT using SIMD
      IntVector v = IntVector.fromArray(SPECIES, x, 0);
      IntVector notV = v.not();

      // Step 2: Add 1 with carry propagation (must be done sequentially)
      notV.intoArray(x, 0);
      int carry = 1;
      for (int i = 0; i < xLen; i++) {
        long temp = (x[i] & MASK_L) + carry;
        x[i] = (int) temp;
        carry = (int) (temp >>> 32);
        if (carry == 0) break;
      }
    } else {
      negate(x, xLen);
    }
  }

  private static void absInplace(final int[] x, final int xLen) {
    if (isNeg(x, xLen)) negate(x, xLen);
  }

  private static void absInto(final int[] dst, final int[] src, final int srcLen) {
    System.arraycopy(src, 0, dst, 0, srcLen);
    absInplace(dst, dst.length);
  }

  private static int numberOfLeadingZeros(final int[] x, final int xLen) {
    int leadingIndex = xLen - 1;
    while ((leadingIndex >= 0) && (x[leadingIndex] == 0)) leadingIndex--;
    return 32 * (xLen - leadingIndex - 1) + Integer.numberOfLeadingZeros(x[leadingIndex]);
  }

  private static void shiftLeftInto(
      final int[] result, final int[] x, final int xLen, final int shift) {
    // Unchecked: result should be initialised with zeroes
    // Unchecked: result length should be at least x.length + limbShift
    int limbShift = shift / N_BITS_PER_LIMB;
    int bitShift = shift % N_BITS_PER_LIMB;
    if (bitShift == 0) {
      System.arraycopy(x, 0, result, limbShift, xLen);
      return;
    }

    int j = limbShift;
    int carry = 0;
    for (int i = 0; i < xLen; ++i, ++j) {
      result[j] = (x[i] << bitShift) | carry;
      carry = x[i] >>> (32 - bitShift);
    }
    if (carry != 0) result[j] = carry; // last carry
  }

  private static void shiftRightInto(
      final int[] result, final int[] x, final int xLen, final int shift) {
    // Unchecked: result length should be at least x.length - limbShift
    int limbShift = shift / 32;
    int bitShift = shift % 32;
    int nLimbs = xLen - limbShift;
    if (nLimbs <= 0) return;

    if (bitShift == 0) {
      System.arraycopy(x, limbShift, result, 0, nLimbs);
      return;
    }

    int carry = 0;
    for (int i = nLimbs - 1 + limbShift, j = nLimbs - 1; j >= 0; i--, j--) {
      int r = (x[i] >>> bitShift) | carry;
      result[j] = r;
      carry = x[i] << (32 - bitShift);
    }
  }

  private static int[] addWithCarry(final int[] x, final int xLen, final int[] y, final int yLen) {
    // Step 1: Add with carry
    int[] a;
    int[] b;
    int aLen;
    int bLen;
    if (xLen < yLen) {
      a = y;
      aLen = yLen;
      b = x;
      bLen = xLen;
    } else {
      a = x;
      aLen = xLen;
      b = y;
      bLen = yLen;
    }
    int[] sum = new int[aLen + 1];
    long carry = 0;
    for (int i = 0; i < bLen; i++) {
      long ai = a[i] & MASK_L;
      long bi = b[i] & MASK_L;
      long s = ai + bi + carry;
      sum[i] = (int) s;
      carry = s >>> 32;
    }
    int icarry = (int) carry;
    for (int i = bLen; i < aLen; i++) {
      sum[i] = a[i] + icarry;
      icarry = (a[i] != 0 && sum[i] == 0) ? 1 : 0;
    }
    sum[aLen] = icarry;
    return sum;
  }

  private static int[] addMul(final int[] a, final int aLen, final int[] b, final int bLen) {
    // Shortest in outer loop, swap if needed
    int[] x;
    int xLen;
    int[] y;
    int yLen;
    if (a.length < b.length) {
      x = b;
      xLen = bLen;
      y = a;
      yLen = aLen;
    } else {
      x = a;
      xLen = aLen;
      y = b;
      yLen = bLen;
    }
    int[] lhs = new int[xLen + yLen + 1];

    // Main algo
    for (int i = 0; i < yLen; i++) {
      long carry = 0;
      long yi = y[i] & MASK_L;

      int k = i;
      for (int j = 0; j < xLen; j++, k++) {
        long prod = yi * (x[j] & MASK_L);
        long sum = (lhs[k] & MASK_L) + prod + carry;
        lhs[k] = (int) sum;
        carry = sum >>> 32;
      }

      // propagate leftover carry
      while (carry != 0 && k < lhs.length) {
        long sum = (lhs[k] & MASK_L) + carry;
        lhs[k] = (int) sum;
        carry = sum >>> 32;
        k++;
      }
    }
    return lhs;
  }

  // ============================================================================
  // SIMD Addition Support
  // ============================================================================
  //
  // This section implements SIMD-accelerated 256-bit unsigned integer addition
  // using the JDK Vector API (jdk.incubator.vector). The implementation treats
  // 8 x 32-bit limbs as 4 x 64-bit lanes for hardware-accelerated parallel
  // addition.
  //
  // Hardware Requirements:
  // - x86-64 with AVX2 (256-bit SIMD) - widely available since 2013
  // - ARM with NEON/SVE (alternative vector instruction sets)
  //
  // Performance Benefits:
  // - 4 additions execute in parallel (SIMD lanes)
  // - Reduced pipeline stalls from data dependencies
  // - Efficient carry propagation via lookup table
  // - JIT compiler generates native SIMD instructions (e.g., vpaddd, vpaddq)
  //
  // Algorithm Overview:
  // 1. Load operands into 256-bit vectors (4 x 64-bit lanes)
  // 2. Parallel addition across all lanes simultaneously
  // 3. Detect per-lane carries using unsigned comparison (result < operand)
  // 4. Detect cascade conditions (lanes with all 1s propagate carries)
  // 5. Calculate cross-lane carry propagation (scalar operations)
  // 6. Apply cascaded carries via lookup table
  // 7. Return result with overflow flag
  //
  // Based on the reference implementation from .NET Runtime:
  // https://github.com/dotnet/runtime/blob/main/src/libraries/System.Runtime.Numerics/src/System/UInt256.cs
  // ============================================================================

  // Lookup table for SIMD carry propagation
  // Each entry represents carry values for 4 lanes based on cascade pattern (16 entries x 32 bytes)
  private static final byte[] BROADCAST_LOOKUP = {
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  };

  /**
   * SIMD-accelerated addition of two UInt256 values using JDK Vector API.
   *
   * <p>This implementation uses the JDK Vector API with 256-bit vectors (4 x 64-bit lanes) to
   * perform hardware-accelerated parallel addition with carry propagation across lanes.
   *
   * <p>Requires: --add-modules jdk.incubator.vector
   *
   * @param a First UInt256 operand
   * @param b Second UInt256 operand
   * @return Result UInt256 containing a + b, and whether overflow occurred
   */
  public static int[] addSIMD(final UInt256 a, final UInt256 b) {
    // Convert to 256-bit vectors with 4 x 64-bit lanes
    long[] aLongs = a.toLongLimbs();
    long[] bLongs = b.toLongLimbs();
    LongVector aVec = LongVector.fromArray(SPECIES_L, aLongs, 0);
    LongVector bVec = LongVector.fromArray(SPECIES_L, bLongs, 0);

    // Lanewise add, carry and cascade
    LongVector resultVec = aVec.add(bVec);
    VectorMask<Long> carryMask = resultVec.compare(VectorOperators.UNSIGNED_LT, aVec);
    VectorMask<Long> cascadeMask = resultVec.compare(VectorOperators.EQ, ALLSET_VEC);

    // Calculate cross-lane carries using scalar operations
    int carry = maskToInt(carryMask);
    int cascade = maskToInt(cascadeMask);
    carry = cascade + 2 * carry;
    cascade ^= carry;
    cascade &= 0x0f;

    int lookupOffset = cascade * 32;
    long[] cascadedCarries = new long[4];
    for (int i = 0; i < 4; i++) {
      long value = 0;
      for (int j = 0; j < 8; j++) {
        value |= ((long) (BROADCAST_LOOKUP[lookupOffset + i * 8 + j] & 0xFF)) << (j * 8);
      }
      cascadedCarries[i] = value;
    }

    // Apply cascaded carries
    LongVector carryVector = LongVector.fromArray(SPECIES_L, cascadedCarries, 0);
    resultVec = resultVec.add(carryVector);
    long[] resultLongs = new long[4];
    resultVec.intoArray(resultLongs, 0);

    int[] result = new int[N_LIMBS + 1];
    for (int i = 0; i < 4; i++) {
      result[i * 2] = (int) resultLongs[i];
      result[i * 2 + 1] = (int) (resultLongs[i] >>> 32);
    }
    boolean overflow = (carry & 0b1_0000) != 0;
    result[N_LIMBS] = overflow ? 1 : 0;
    return result;
  }

  /**
   * Convert a VectorMask to an integer bitmask.
   *
   * <p>Each bit represents whether the corresponding lane's mask is true.
   *
   * @param mask The vector mask to convert
   * @return Integer bitmask where bit i is set if lane i's mask is true
   */
  private static int maskToInt(final VectorMask<Long> mask) {
    int result = 0;
    for (int i = 0; i < 4; i++) {
      if (mask.laneIsSet(i)) {
        result |= (1 << i);
      }
    }
    return result;
  }

  /**
   * SIMD-accelerated addition using IntVector (8 x 32-bit lanes).
   *
   * <p>This implementation uses IntVector for a more direct approach, performing parallel addition
   * across all 8 limbs, then handling carry propagation sequentially.
   *
   * <p>Simpler than the LongVector approach but requires sequential carry propagation.
   *
   * @param a First UInt256 operand
   * @param b Second UInt256 operand
   * @return Result UInt256 containing a + b, and whether overflow occurred
   */
  public static int[] addSIMDInt(final UInt256 a, final UInt256 b) {
    IntVector aVec = IntVector.fromArray(SPECIES, a.limbs, 0);
    IntVector bVec = IntVector.fromArray(SPECIES, b.limbs, 0);
    IntVector sumVec = aVec.add(bVec);
    VectorMask<Integer> carryMask = sumVec.compare(VectorOperators.UNSIGNED_LT, aVec);

    int[] result = new int[N_LIMBS + 1];
    sumVec.intoArray(result, 0);
    long carry = 0;
    for (int i = 0; i < N_LIMBS; i++) {
      if (carryMask.laneIsSet(i)) {
        carry = 1;
      }

      if (i < N_LIMBS - 1 && carry != 0) {
        long next = (result[i + 1] & MASK_L) + carry;
        result[i + 1] = (int) next;
        carry = next >>> 32;
      }
    }
    boolean overflow = carry != 0;
    result[N_LIMBS] = overflow ? 1 : 0;
    return result;
  }

  private static int[] knuthRemainder(final int[] dividend, final int[] modulus) {
    int[] result = new int[N_LIMBS];
    int divLen = nSetLimbs(dividend);
    int modLen = nSetLimbs(modulus);
    int cmp = compareLimbs(dividend, divLen, modulus, modLen);
    if (cmp < 0) {
      System.arraycopy(dividend, 0, result, 0, divLen);
      return result;
    } else if (cmp == 0) {
      return result;
    }

    int shift = numberOfLeadingZeros(modulus, modLen);
    int limbShift = shift / 32;
    int n = modLen - limbShift;
    if (n == 0) return result;
    if (n == 1) {
      if (divLen == 1) {
        result[0] = Integer.remainderUnsigned(dividend[0], modulus[0]);
        return result;
      }
      long d = modulus[0] & MASK_L;
      long rem = 0;
      // Process from most significant limb downwards
      for (int i = divLen - 1; i >= 0; i--) {
        long cur = (rem << 32) | (dividend[i] & MASK_L);
        rem = Long.remainderUnsigned(cur, d);
      }
      result[0] = (int) rem;
      result[1] = (int) (rem >>> 32);
      return result;
    }
    // Normalize
    int m = divLen - n;
    int bitShift = shift % 32;
    int[] vLimbs = new int[n];
    shiftLeftInto(vLimbs, modulus, modLen, shift);
    int[] uLimbs = new int[divLen + 1];
    shiftLeftInto(uLimbs, dividend, divLen, bitShift);

    long[] vLimbsAsLong = new long[n];
    for (int i = 0; i < n; i++) {
      vLimbsAsLong[i] = vLimbs[i] & MASK_L;
    }

    // Main division loop
    long vn1 = vLimbsAsLong[n - 1];
    long vn2 = vLimbsAsLong[n - 2];
    for (int j = m; j >= 0; j--) {
      long ujn = (uLimbs[j + n] & MASK_L);
      long ujn1 = (uLimbs[j + n - 1] & MASK_L);
      long ujn2 = (uLimbs[j + n - 2] & MASK_L);

      long dividendPart = (ujn << 32) | ujn1;
      // Check that no need for Unsigned version of divrem.
      long qhat = Long.divideUnsigned(dividendPart, vn1);
      long rhat = Long.remainderUnsigned(dividendPart, vn1);

      while (qhat == 0x1_0000_0000L || Long.compareUnsigned(qhat * vn2, (rhat << 32) | ujn2) > 0) {
        qhat--;
        rhat += vn1;
        if (rhat >= 0x1_0000_0000L) break;
      }

      // Multiply-subtract qhat*v from u slice
      long borrow = 0;
      for (int i = 0; i < n; i++) {
        long prod = vLimbsAsLong[i] * qhat;
        long sub = (uLimbs[i + j] & MASK_L) - (prod & MASK_L) - borrow;
        uLimbs[i + j] = (int) sub;
        borrow = (prod >>> 32) - (sub >> 32);
      }
      long sub = (uLimbs[j + n] & MASK_L) - borrow;
      uLimbs[j + n] = (int) sub;

      if (sub < 0) {
        // Add back
        long carry = 0;
        for (int i = 0; i < n; i++) {
          long sum = (uLimbs[i + j] & MASK_L) + vLimbsAsLong[i] + carry;
          uLimbs[i + j] = (int) sum;
          carry = sum >>> 32;
        }
        uLimbs[j + n] = (int) (uLimbs[j + n] + carry);
      }
    }
    // Unnormalize remainder
    shiftRightInto(result, uLimbs, n, bitShift);
    return result;
  }
  // --------------------------------------------------------------------------
  // endregion
}
