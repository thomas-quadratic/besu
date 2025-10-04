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

/**
 * 256-bits wide (unsigned) integer class.
 *
 * <p>This class is an optimised version of BigInteger for fixed width 256-bits integers.
 * In contrast to Java int, by default operations interpretes UInt256 as unsigned, with specific signed operations when needed.
 */
public final class UInt256 {
  // region Internals
  // --------------------------------------------------------------------------
  // Internally, it is represented as 4 long limbs in little-endian order, that is from least to most significant limbs.
  // Note that Java's long are themselves big-endian.
  //
  // Length is used to optimise algorithms, skipping leading zeroes. It should be by default 4, but can be less in cases we know for sure that the 4 - length most significant limbs are zero. This gives optional flexibility for optimisations, for example in loops. In all cases, all 4 limbs are allocated and length can be safely ignored.

  /** Fixed size in bytes. */
  public static final int BYTESIZE = 32;

  /** Fixed size in bits. */
  public static final int BITSIZE = 256;

  // Fixed number of limbs or digits
  private static final int N_LIMBS = 4;
  // Fixed number of bytes per limb.
  private static final int N_BYTES_PER_LIMB = 8;
  // Fixed number of bits per limb.
  private static final int N_BITS_PER_LIMB = 64;
  // Mask for unsigned int -> long conversion
  private static final long MASK_L = 0xFFFFFFFFL;

  private final long[] limbs;
  private final int length;

  long[] limbs() {
    return limbs;
  }

  // --------------------------------------------------------------------------
  // endregion

  /** The constant 0. */
  public static final UInt256 ZERO = new UInt256(new long[] {0, 0, 0, 0}, 0);

  // region Constructors
  // --------------------------------------------------------------------------

  UInt256(final long[] limbs, final int length) {
    // Unchecked length: assumes limbs have length == N_LIMBS
    this.limbs = limbs;
    this.length = length;
  }

  UInt256(final long[] limbs) {
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
    long[] limbs = new long[N_LIMBS];
    for (int i = nLimbs - 1; i >= 0; i--) {
      limbs[i] = buf.getLong();
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
    long[] limbs = new long[N_LIMBS];
    limbs[0] = value & MASK_L;
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
    long[] limbs = new long[N_LIMBS];
    limbs[0] = value;
    return new UInt256(limbs, 1);
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
  public static UInt256 fromArray(final long[] arr) {
    long[] limbs = new long[N_LIMBS];
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
    return (int) limbs[0];
  }

  /**
   * Convert to long.
   *
   * @return Value truncated to a long, possibly lossy.
   */
  public long longValue() {
    return limbs[0];
  }

  /**
   * Convert to BigEndian byte array.
   *
   * @return Big-endian ordered bytes for this UInt256 value.
   */
  public byte[] toBytesBE() {
    ByteBuffer buf = ByteBuffer.allocate(BYTESIZE).order(ByteOrder.BIG_ENDIAN);
    for (int i = N_LIMBS - 1; i >= 0; i--) {
      buf.putLong(limbs[i]);
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
    return (limbs[0] | limbs[1] | limbs[2] | limbs[3]) == 0;
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
      comp = Long.compareUnsigned(a.limbs[i], b.limbs[i]);
      if (comp != 0) return comp;
    }
    return 0;
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) return true;
    if (!(obj instanceof UInt256)) return false;
    UInt256 other = (UInt256) obj;

    long xor =
        (this.limbs[0] ^ other.limbs[0])
            | (this.limbs[1] ^ other.limbs[1])
            | (this.limbs[2] ^ other.limbs[2])
            | (this.limbs[3] ^ other.limbs[3]);
    return xor == 0;
  }

  @Override
  public int hashCode() {
    int h = 1;
    for (int i = 0; i < N_LIMBS; i++) {
      h = 31 * h + Long.hashCode(limbs[i]);
    }
    return h;
  }

  /**
   * Is the two complements signed representation of this integer negative.
   *
   * @return True if the two complements representation of this integer is negative.
   */
  public boolean isNegative() {
    return (length == 4) && (isNeg(limbs, 4));
  }

  // --------------------------------------------------------------------------
  // endregion

  // region Bitwise Operations
  // --------------------------------------------------------------------------

  /**
   * Shifts value to the left.
   *
   * @param shift number of bits to shift. If negative, shift right instead.
   * @return Shifted UInt256 value.
   */
  public UInt256 shiftLeft(final int shift) {
    if (shift < 0) return shiftRight(-shift);
    if (shift >= BITSIZE) return ZERO;
    if (shift == 0 || isZero()) return this;
    long[] shifted = new long[N_LIMBS];
    shiftLeftInto(shifted, this.limbs, this.length, shift);
    return new UInt256(shifted);
  }

  /**
   * Shifts value to the right.
   *
   * @param shift number of bits to shift. If negative, shift left instead.
   * @return Shifted UInt256 value.
   */
  public UInt256 shiftRight(final int shift) {
    if (shift < 0) return shiftLeft(-shift);
    if (shift >= length * N_BITS_PER_LIMB) return ZERO;
    if (shift == 0 || isZero()) return this;
    long[] shifted = new long[N_LIMBS];
    shiftRightInto(shifted, this.limbs, this.length, shift);
    return new UInt256(shifted);
  }

  // --------------------------------------------------------------------------
  // endregion

  // region Arithmetic Operations
  // --------------------------------------------------------------------------

  /**
   * (Signed) absolute value
   *
   * @return The absolute value of this signed integer.
   */
  public UInt256 abs() {
    long[] newLimbs = Arrays.copyOf(limbs, N_LIMBS);
    absInplace(newLimbs, N_LIMBS);
    return new UInt256(newLimbs);
  }

  /**
   * Addition (modulo 2**256).
   *
   * @param other The integer to add.
   * @return The sum of this with other.
   */
  public UInt256 add(final UInt256 other) {
    if (other.isZero()) return this;
    if (this.isZero()) return other;
    long[] newLimbs = addWithCarry(this.limbs, this.length, other.limbs, other.length);
    return new UInt256(newLimbs, Math.min(N_LIMBS, newLimbs.length));
  }

  /**
   * Multiplication (modulo 2**256).
   *
   * @param other The integer to add.
   * @return The sum of this with other.
   */
  public UInt256 mul(final UInt256 other) {
    if (this.isZero() || other.isZero()) return ZERO;
    long[] prod = addMul(this.limbs, this.length, other.limbs, other.length);
    long[] result = Arrays.copyOf(prod, N_LIMBS);
    return new UInt256(result, Math.min(N_LIMBS, result.length));
  }

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
   * @param modulus The modulus of the reduction
   * @return The remainder modulo {@code modulus}.
   */
  public UInt256 signedMod(final UInt256 modulus) {
    if (this.isZero() || modulus.isZero()) return ZERO;
    long[] x = new long[N_LIMBS];
    long[] y = new long[N_LIMBS];
    absInto(x, this.limbs, N_LIMBS);
    absInto(y, modulus.limbs, N_LIMBS);
    long[] r = knuthRemainder(x, y);
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
    // if (this.isZero()) return other.mod(modulus);
    // if (other.isZero()) return this.mod(modulus);
    long[] sum = addWithCarry(this.limbs, this.length, other.limbs, other.length);
    long[] rem = knuthRemainder(sum, modulus.limbs);
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
    long[] result = addMul(this.limbs, this.length, other.limbs, other.length);
    result = knuthRemainder(result, modulus.limbs);
    return new UInt256(result, modulus.length);
  }

  // --------------------------------------------------------------------------
  // endregion

  // region Support (private) Algorithms
  // --------------------------------------------------------------------------
  private static int nLeadingZeroBits(final long[] x, final int xLen) {
    int leadingIndex = xLen - 1;
    while ((leadingIndex >= 0) && (x[leadingIndex] == 0)) leadingIndex--;
    return N_BITS_PER_LIMB * (xLen - leadingIndex - 1) + Long.numberOfLeadingZeros(x[leadingIndex]);
  }

  private static int nLeadingZeroLimbs(final long[] x, final int maxLength) {
    int offset = maxLength - 1;
    while ((offset >= 0) && (x[offset] == 0)) offset--;
    return offset + 1;
  }

  private static int compareLimbs(final long[] a, final int aLen, final long[] b, final int bLen) {
    int cmp;
    if (aLen > bLen) {
      for (int i = aLen - 1; i >= bLen; i--) {
        if (a[i] != 0) return 1;
      }
    } else if (aLen < bLen) {
      for (int i = bLen - 1; i >= aLen; i--) {
        if (b[i] != 0) return -1;
      }
    }
    for (int i = Math.min(aLen, bLen) - 1; i >= 0; i--) {
      cmp = Long.compareUnsigned(a[i], b[i]);
      if (cmp != 0) return cmp;
    }
    return 0;
  }

  private static boolean isNeg(final long[] x, final int xLen) {
    return x[xLen - 1] < 0;
  }

  private static void negate(final long[] x, final int xLen) {
    int carry = 1;
    for (int i = 0; i < xLen; i++) {
      x[i] = ~x[i] + carry;
      carry = (x[i] == 0 && carry == 1) ? 1 : 0;
    }
  }

  private static void absInplace(final long[] x, final int xLen) {
    if (isNeg(x, xLen)) negate(x, xLen);
  }

  private static void absInto(final long[] dst, final long[] src, final int srcLen) {
    System.arraycopy(src, 0, dst, 0, srcLen);
    absInplace(dst, dst.length);
  }

  private static void shiftLeftInto(
      final long[] result, final long[] x, final int xLen, final int shift) {
    // Unchecked: result should be initialised with zeroes
    // Unchecked: result length should be at least x.length + limbShift + 0 or 1 for carry
    int limbShift = shift / N_BITS_PER_LIMB;
    int bitShift = shift % N_BITS_PER_LIMB;
    if (limbShift >= xLen) return;
    if (bitShift == 0) {
      System.arraycopy(x, 0, result, limbShift, xLen);
      return;
    }

    int j = limbShift;
    long carry = 0;
    for (int i = 0; i < xLen; ++i, ++j) {
      result[j] = (x[i] << bitShift) | carry;
      carry = x[i] >>> (N_BITS_PER_LIMB - bitShift);
    }
    if (carry != 0) result[j] = carry; // last carry
  }

  private static void shiftRightInto(
      final long[] result, final long[] x, final int xLen, final int shift) {
    // Unchecked: result length should be at least x.length - limbShift
    int limbShift = shift / N_BITS_PER_LIMB;
    int bitShift = shift % N_BITS_PER_LIMB;
    int nLimbs = xLen - limbShift;
    if (nLimbs <= 0) return;

    if (bitShift == 0) {
      System.arraycopy(x, limbShift, result, 0, nLimbs);
      return;
    }

    long carry = 0;
    for (int i = nLimbs - 1 + limbShift, j = nLimbs - 1; j >= 0; i--, j--) {
      long r = (x[i] >>> bitShift) | carry;
      result[j] = r;
      carry = x[i] << (N_BITS_PER_LIMB - bitShift);
    }
  }

  private static long[] addWithCarry(final long[] x, final int xLen, final long[] y, final int yLen) {
    // Step 1: Add with carry
    long[] a;
    long[] b;
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
    long[] sum = new long[aLen + 1];
    long carry = 0;
    for (int i = 0; i < bLen; i++) {
      long s = a[i] + carry;
      carry = (Long.compareUnsigned(s, a[i]) < 0) ? 1 : 0;
      s += b[i];
      if (Long.compareUnsigned(s, b[i]) < 0) carry++;
    }
    int i = bLen;
    while ((carry != 0) && (i < aLen + 1)) {
      sum[i] = a[i] + carry;
      carry = (Long.compareUnsigned(sum[i], carry) < 0) ? 1 : 0;
      i++;
    }
    for (; i < aLen; i++) sum[i] = a[i];
    return sum;
  }

  private static long[] addMul(final long[] a, final int aLen, final long[] b, final int bLen) {
    // Shortest in outer loop, swap if needed
    long[] x;
    int xLen;
    long[] y;
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
    long[] lhs = new long[xLen + yLen + 1];
    for (int i = 0; i < y.length; i++) {
      long carry = 0;
      int k = i;
      for (int j = 0; j < x.length; j++, k++) {
        long p0 = y[i] * x[j];
        long p1 = Math.multiplyHigh(y[i], x[j]);
        p0 += carry;
        if (Long.compareUnsigned(p0, carry) < 0) p1++;
        lhs[k] += p0;
        if (Long.compareUnsigned(lhs[k], p0) < 0) p1++;
        carry = p1;
      }

      // propagate leftover carry
      while (carry != 0 && k < lhs.length) {
        lhs[k] += carry;
        carry = (Long.compareUnsigned(lhs[k], carry) < 0) ? 1 : 0;
        k++;
      }
    }
    return lhs;
  }

  private static long reciprocal(final long x) {
    return x;
  }

  private static long reciprocal(final long x0, final long x1) {
    return x0 + x1;
  }

  private static long div2by1(final long[] x, final int xLen, final long y, final long yInv) {
    // wrapping umul x1 * yInv
    long q0 = x[xLen - 1] * yInv;
    long q1 = Math.multiplyHigh(x[xLen - 1], yInv);

    // wrapping uadd <q1, q0> + <x1, x0> + <1, 0>
    q0 += x[xLen - 2];
    long carry = (Long.compareUnsigned(q0, x[xLen - 2]) < 0) ? 1 : 0;
    q1 += x[xLen -1] + carry + 1;

    x[xLen -2] -= q1 * y;

    if (Long.compareUnsigned(x[xLen - 2], q0) > 0) {
      q1 -= 1;
      x[xLen - 2] += y;
    }

    if (Long.compareUnsigned(x[xLen - 2], y) >= 0) {
      q1 += 1;
      x[xLen - 2] -= y;
    }
    x[xLen - 1] = 0;
    return q1;
  }

  private static long div3by2(final long[] x, final int xLen, final long y0, final long y1, final long yInv) {
    // <x2, x1, x0> divided by <y1, y0>.
    // Works inplace on x, modifying it to hold the remainder.
    // Returns the quotient q.
    // Requires <x2, x1> < <y1, y0> otherwise quotient overflows.
    long overflow;  // carry or borrow
    long x0 = x[xLen - 3];
    long x1 = x[xLen - 2];
    long x2 = x[xLen - 1];

    // wrapping umul x2 * yInv
    long q0 = x2 * yInv;
    long q1 = Math.multiplyHigh(x2, yInv);

    // wrapping uadd <q1, q0> + <x2, x1>
    q0 = q0 + x1;
    overflow = (Long.compareUnsigned(q0, x1) < 0) ? 1 : 0;
    q1 = q1 + x2 + overflow;
    
    x[xLen - 2] -= q1 * y1;

    // wrapping umul q1 * y0
    long t0 = q1 * y0;
    long t1 = Math.multiplyHigh(q1, y0);

    // wrapping sub <r1, x0> − <t1, t0> − <y1, y0>
    overflow = Long.compareUnsigned(x0, t0) < 0 ? 1 : 0;
    x[xLen - 3] -= t0;
    x[xLen - 2] -= (t1 + overflow);

    overflow = Long.compareUnsigned(x[xLen - 3], y0) < 0 ? 1 : 0;
    x[xLen - 3] -= y0;
    x[xLen - 2] -= (y1 + overflow);

    q1 += 1;
    if (x[xLen - 2] >= q0) {
        q1 -= 1;
	x[xLen - 3] += y0;
        overflow = (Long.compareUnsigned(x[xLen - 3], y0) < 0) ? 1 : 0;
	x[xLen - 2] += y1 + overflow;
    }

    int cmp = Long.compareUnsigned(x[xLen - 2], y1);
    if ((cmp > 0) || ((cmp == 0) && (Long.compareUnsigned(x[xLen - 3], y0) >= 0))) {
      q1 += 1;
      overflow = Long.compareUnsigned(x[xLen - 3], y0) < 0 ? 1 : 0;
      x[xLen - 3] -= y0;
      x[xLen - 2] -= (y1 + overflow);
    }
    return q1;
  }

  private static void modNby1(final long[] x, final int xLen, final long y, final long yInv) {
    for (int i = xLen - 1; i >= 1; i--) {
      // if x >= y, overflows, so must substract y first
      if (Long.compareUnsigned(x[i], y) >= 0) x[i] -= y;
      div2by1(x, i + 1, y, yInv);
    }
  }

  private static void modNby2(final long[] x, final int xLen, final long y0, final long y1, final long yInv) {
    for (int i = xLen - 1; i >= 1; i--) {
      // if <x+1, x+0> >= <y1, y0>, overflows, so must substract <y1, y0> first
      int cmp = Long.compareUnsigned(x[i], y1);
      if ((cmp > 0) || ((cmp == 0) && (Long.compareUnsigned(x[i - 1], y0) >= 0))) {
	long borrow = (Long.compareUnsigned(y0, x[i - 1]) > 0) ? 1 : 0;
	x[i - 1] -= y0;
        x[i] -= (y1 + borrow);
      }
      div3by2(x, i + 1, y0, y1, yInv);
    }
  }

  private static void modNbyM(final long[] x, final int xLen, final long[] y, final int yLen, final long yInv) {
    long yn1 = y[yLen - 1];
    long yn0 = y[yLen - 2];
    int m = xLen - yLen;
    for (int j = m - 1; j >= 0; j--) {
      long borrow = 0;

      // Unlikely overflow 3x2 case: x[..2] == y
      if (((yn1 ^ x[j + yLen]) | (yn0 ^ x[j + yLen - 1])) == 0) {
        // q = <0, 1> in this case, so multiply is trivial.
	// Still need to substract (and add back if needed).
        for (int i = 0; i < yLen; i++) {
          x[j + i + 1] -= borrow;
          borrow = (Long.compareUnsigned(x[i + j], borrow) < 0) ? 1 : 0;
          x[j + i + 1] -= y[i];
          if (Long.compareUnsigned(x[j + i + 1], y[i]) < 0) borrow++;
        }
        x[j + yLen] -= borrow;
        borrow = (Long.compareUnsigned(x[j + yLen], borrow) < 0) ? 1 : 0;
      } else {
        long q = div3by2(x, j + yLen + 1, yn0, yn1, yInv);

        // Multiply-subtract: already have highest 2 limbs
        for (int i = 0; i < yLen - 2; i++) {
          long p0 = y[i] * q;
          long p1 = Math.multiplyHigh(y[i], q);
          x[j + i] -= borrow;
          borrow = (Long.compareUnsigned(x[j + i], borrow) < 0) ? p1 + 1 : p1;
          x[j + i] -= p0;
          if (Long.compareUnsigned(x[j + i], p0) < 0) borrow++;
        }
        x[j + yLen - 1] -= borrow;
        borrow = (Long.compareUnsigned(x[j + yLen - 1], borrow) < 0) ? 1 : 0;
        x[j + yLen] -= borrow;
        borrow = (Long.compareUnsigned(x[j + yLen], borrow) < 0) ? 1 : 0;
      }

      if (borrow != 0) { // unlikely
        // Add back
        long carry = 0;
	int i;
        for (i = 0; i < yLen; i++) {
          x[j + i] = x[j + i] + carry;
          carry = (Long.compareUnsigned(x[j + i], carry) < 0) ? 1 : 0;
          x[j + i] = x[j + i] + y[i];
          if (Long.compareUnsigned(x[j + i], carry) < 0) carry++;
        }
	while ((carry != 0) && (i < xLen)) {
          x[j + i] = x[j + i] + carry;
          carry = (Long.compareUnsigned(x[j + i], carry) < 0) ? 1 : 0;
        }
      }
    }
  }

  private static long[] knuthRemainder(final long[] dividend, final long[] modulus) {
    long[] result = new long[N_LIMBS];
    int divLen = dividend.length - nLeadingZeroLimbs(dividend, dividend.length);
    int modLen = modulus.length - nLeadingZeroLimbs(modulus, modulus.length);

    // Shortcuts
    if (modLen == 0) return result;
    if (divLen == 1) {
      if (modLen == 1) result[0] = Long.remainderUnsigned(dividend[0], modulus[0]);
      else System.arraycopy(dividend, 0, result, 0, divLen);
      return result;
    }
    int cmp = compareLimbs(dividend, divLen, modulus, modLen);
    if (cmp < 0) {
      System.arraycopy(dividend, 0, result, 0, divLen);
      return result;
    } else if (cmp == 0) {
      return result;
    }

    // Perform Division
    // -- Normalise
    int shift = nLeadingZeroBits(modulus, modLen);
    int uLen = divLen + 1;
    long[] uLimbs = new long[uLen];
    shiftLeftInto(uLimbs, dividend, divLen, shift);
    long[] vLimbs = new long[modLen];
    shiftLeftInto(vLimbs, modulus, modLen, shift);

    // -- Divide
    if (modLen == 1) {
      long inv = reciprocal(vLimbs[0]);
      modNby1(uLimbs, uLen, vLimbs[0], inv);
    } else if (modLen == 2) {
      long inv = reciprocal(vLimbs[0], vLimbs[1]);
      modNby2(uLimbs, uLen, vLimbs[0], vLimbs[1], inv);
    } else {
      long inv = reciprocal(vLimbs[modLen - 2], vLimbs[modLen - 1]);
      modNbyM(uLimbs, uLen, vLimbs, modLen, inv);
    }
    // -- Unnormalize
    shiftRightInto(result, uLimbs, modLen, shift);
    return result;
  }
    
  // --------------------------------------------------------------------------
  // endregion
}
