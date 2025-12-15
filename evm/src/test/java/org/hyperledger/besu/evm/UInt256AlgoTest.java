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

import static org.assertj.core.api.Assertions.assertThat;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;

import java.math.BigInteger;
import java.util.Arrays;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

public class UInt256AlgoTest {

  private byte[] bigIntToBytes(final BigInteger x) {
    byte[] a = x.toByteArray();
    if (a.length == 0) return new byte[0];
    // Remove leading zero byte added by BigInteger for sign
    if (a.length > 32 && a[0] == 0) {
      byte[] result = new byte[32];
      System.arraycopy(a, a.length - 32, result, 0, 32);
      return result;
    }
    // Remove leading zeros
    int firstNonZero = 0;
    while (firstNonZero < a.length && a[firstNonZero] == 0) {
      firstNonZero++;
    }
    if (firstNonZero == a.length) return new byte[0];
    if (firstNonZero == 0) return a;
    byte[] result = new byte[a.length - firstNonZero];
    System.arraycopy(a, firstNonZero, result, 0, result.length);
    return result;
  }

  private byte[] bigIntToBytes32(final BigInteger x) {
    BigInteger mod = x.mod(BigInteger.ONE.shiftLeft(256));
    byte[] bytes = mod.toByteArray();
    byte[] result = new byte[32];
    if (bytes.length == 0) return result;
    if (bytes.length <= 32) {
      System.arraycopy(bytes, 0, result, 32 - bytes.length, bytes.length);
    } else {
      System.arraycopy(bytes, bytes.length - 32, result, 0, 32);
    }
    return result;
  }

  // region Addition Tests
  // Note: addIntWidening, addIntAndCarry convert variable-length input to 8 int limbs internally, but accept any length
  // addByteVarLen works directly with variable-length byte arrays
  // addByteFixedLen expects exactly 32-byte inputs

  @Test
  public void testAddA_simple() {
    byte[] a = new byte[] {5};
    byte[] b = new byte[] {10};
    byte[] result = UInt256Algo.addIntWidening(a, b);
    byte[] expected = new byte[] {15};
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testAddA_withCarry() {
    byte[] a = new byte[] {(byte) 0xFF};
    byte[] b = new byte[] {2};
    byte[] result = UInt256Algo.addIntWidening(a, b);
    byte[] expected = new byte[] {1, 1}; // 255 + 2 = 257 = 0x0101
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testAddA_zeroLeft() {
    byte[] a = new byte[0];
    byte[] b = new byte[] {42};
    byte[] result = UInt256Algo.addIntWidening(a, b);
    assertArrayEquals(b, result);
  }

  @Test
  public void testAddA_zeroRight() {
    byte[] a = new byte[] {42};
    byte[] b = new byte[0];
    byte[] result = UInt256Algo.addIntWidening(a, b);
    assertArrayEquals(a, result);
  }

  @Test
  public void testAddA_overflow() {
    byte[] maxValue = new byte[32];
    Arrays.fill(maxValue, (byte) 0xFF);
    byte[] one = new byte[] {1};
    byte[] result = UInt256Algo.addIntWidening(maxValue, one);
    // Should wrap to zero (modulo 2^256)
    byte[] expected = new byte[0];
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testAddB_simple() {
    byte[] a = new byte[] {7};
    byte[] b = new byte[] {3};
    byte[] result = UInt256Algo.addIntAndCarry(a, b);
    byte[] expected = new byte[] {10};
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testAddC_simple() {
    byte[] a = new byte[] {20};
    byte[] b = new byte[] {30};
    byte[] result = UInt256Algo.addByteVarLen(a, b);
    byte[] expected = new byte[] {50};
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testAddC_withCarry() {
    byte[] a = new byte[] {(byte) 0xFF, (byte) 0xFF};
    byte[] b = new byte[] {1};
    byte[] result = UInt256Algo.addByteVarLen(a, b);
    byte[] expected = new byte[] {1, 0, 0}; // 0xFFFF + 1 = 0x10000
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testAddD_simple() {
    // addByteFixedLen returns exactly 32 bytes always
    byte[] a = bigIntToBytes32(BigInteger.valueOf(100));
    byte[] b = bigIntToBytes32(BigInteger.valueOf(50));
    byte[] result = UInt256Algo.addByteFixedLen(a, b);
    byte[] expected = bigIntToBytes32(BigInteger.valueOf(150));
    assertArrayEquals(expected, result);
  }

  @Test
  public void testAddD_overflow() {
    // addByteFixedLen returns exactly 32 bytes, wrapping on overflow
    byte[] maxValue = new byte[32];
    Arrays.fill(maxValue, (byte) 0xFF);
    byte[] one = bigIntToBytes32(BigInteger.ONE);
    byte[] result = UInt256Algo.addByteFixedLen(maxValue, one);
    byte[] expected = new byte[32]; // All zeros
    assertArrayEquals(expected, result);
  }

  @Test
  public void testAdd_implementations_consistency() {
    // Test that addIntWidening, addIntAndCarry, addByteVarLen produce equivalent results
    byte[] a = new byte[] {(byte) 0x12, (byte) 0x34, (byte) 0x56};
    byte[] b = new byte[] {(byte) 0xAB, (byte) 0xCD};

    byte[] resultA = UInt256Algo.addIntWidening(a, b);
    byte[] resultB = UInt256Algo.addIntAndCarry(a, b);
    byte[] resultC = UInt256Algo.addByteVarLen(a, b);

    assertThat(UInt256Algo.compare(resultA, resultB)).isEqualTo(0);
    assertThat(UInt256Algo.compare(resultA, resultC)).isEqualTo(0);
  }

  @Test
  public void testAdd_largeNumbers() {
    // Test with larger byte arrays
    BigInteger aBig = new BigInteger("123456789abcdef0123456789abcdef", 16);
    BigInteger bBig = new BigInteger("fedcba9876543210fedcba987654321", 16);
    byte[] a = bigIntToBytes(aBig);
    byte[] b = bigIntToBytes(bBig);

    byte[] resultA = UInt256Algo.addIntWidening(a, b);
    byte[] resultB = UInt256Algo.addIntAndCarry(a, b);
    byte[] resultC = UInt256Algo.addByteVarLen(a, b);

    byte[] expected = bigIntToBytes(aBig.add(bBig).mod(BigInteger.ONE.shiftLeft(256)));

    assertThat(UInt256Algo.compare(resultA, expected)).isEqualTo(0);
    assertThat(UInt256Algo.compare(resultB, expected)).isEqualTo(0);
    assertThat(UInt256Algo.compare(resultC, expected)).isEqualTo(0);
  }

  // endregion

  // region Multiplication Tests

  @Test
  public void testMul_simple() {
    byte[] a = new byte[] {5};
    byte[] b = new byte[] {7};
    byte[] result = UInt256Algo.mul(a, b);
    byte[] expected = new byte[] {35};
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testMul_zero() {
    byte[] a = new byte[] {42};
    byte[] zero = new byte[0];
    byte[] result = UInt256Algo.mul(a, zero);
    assertArrayEquals(new byte[0], result);
  }

  @Test
  public void testMul_one() {
    byte[] a = new byte[] {1, 2, 3, 4};
    byte[] one = new byte[] {1};
    byte[] result = UInt256Algo.mul(a, one);
    assertThat(UInt256Algo.compare(result, a)).isEqualTo(0);
  }

  @Test
  public void testMul_large() {
    BigInteger aBig = new BigInteger("123456789abcdef", 16);
    BigInteger bBig = new BigInteger("fedcba987654321", 16);
    byte[] a = bigIntToBytes(aBig);
    byte[] b = bigIntToBytes(bBig);
    byte[] result = UInt256Algo.mul(a, b);
    byte[] expected = bigIntToBytes(aBig.multiply(bBig).mod(BigInteger.ONE.shiftLeft(256)));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  // endregion

  // region Modulo Tests

  @Test
  public void testMod_simple() {
    byte[] a = new byte[] {23};
    byte[] b = new byte[] {5};
    byte[] result = UInt256Algo.mod(a, b);
    byte[] expected = new byte[] {3}; // 23 % 5 = 3
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testMod_dividendSmallerThanModulus() {
    byte[] a = new byte[] {5};
    byte[] b = new byte[] {10};
    byte[] result = UInt256Algo.mod(a, b);
    byte[] expected = new byte[] {5};
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testMod_dividendEqualsModulus() {
    byte[] a = new byte[] {42};
    byte[] b = new byte[] {42};
    byte[] result = UInt256Algo.mod(a, b);
    assertArrayEquals(new byte[0], result);
  }

  @Test
  public void testMod_byZero() {
    byte[] a = new byte[] {42};
    byte[] zero = new byte[0];
    byte[] result = UInt256Algo.mod(a, zero);
    assertArrayEquals(new byte[0], result);
  }

  @Test
  public void testMod_zeroModAnything() {
    byte[] zero = new byte[0];
    byte[] b = new byte[] {42};
    byte[] result = UInt256Algo.mod(zero, b);
    assertArrayEquals(new byte[0], result);
  }

  @Test
  public void testMod_largeNumbers() {
    BigInteger aBig = new BigInteger("cea0c5cc171fa61277e5604a3bc8aef4de3d3882", 16);
    BigInteger bBig = new BigInteger("7dae7454bb193b1c28e64a6a935bc3", 16);
    byte[] a = bigIntToBytes(aBig);
    byte[] b = bigIntToBytes(bBig);
    byte[] result = UInt256Algo.mod(a, b);
    byte[] expected = bigIntToBytes(aBig.mod(bBig));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testMod_singleLimb() {
    byte[] a = new byte[] {1};
    byte[] b = new byte[] {(byte) 0xFF};
    byte[] result = UInt256Algo.mod(a, b);
    byte[] expected = new byte[] {1};
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  // endregion

  // region Signed Modulo Tests

  @Test
  public void testSignedMod_positiveNumbers() {
    byte[] a = new byte[32];
    byte[] b = new byte[32];
    a[31] = 23; // Small positive
    b[31] = 5;
    byte[] result = UInt256Algo.signedMod(a, b);
    byte[] expected = new byte[32];
    expected[31] = 3; // 23 % 5 = 3
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testSignedMod_byZero() {
    byte[] a = new byte[] {42};
    byte[] zero = new byte[0];
    byte[] result = UInt256Algo.signedMod(a, zero);
    assertArrayEquals(new byte[0], result);
  }

  // endregion

  // region Modular Addition Tests

  @Test
  public void testAddMod_simple() {
    byte[] a = new byte[] {10};
    byte[] b = new byte[] {15};
    byte[] m = new byte[] {7};
    byte[] result = UInt256Algo.addMod(a, b, m);
    byte[] expected = new byte[] {4}; // (10 + 15) % 7 = 25 % 7 = 4
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testAddMod_withZeroModulus() {
    byte[] a = new byte[] {10};
    byte[] b = new byte[] {15};
    byte[] zero = new byte[0];
    byte[] result = UInt256Algo.addMod(a, b, zero);
    assertArrayEquals(new byte[0], result);
  }

  @Test
  public void testAddMod_withZeroAddend() {
    byte[] a = new byte[] {10};
    byte[] zero = new byte[0];
    byte[] m = new byte[] {7};
    byte[] result = UInt256Algo.addMod(a, zero, m);
    byte[] expected = new byte[] {3}; // 10 % 7 = 3
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testAddMod_largeNumbers() {
    BigInteger aBig = new BigInteger("1000000000000000000000000", 16);
    BigInteger bBig = new BigInteger("c350", 16);
    BigInteger mBig = new BigInteger("3e8", 16);
    byte[] a = bigIntToBytes(aBig);
    byte[] b = bigIntToBytes(bBig);
    byte[] m = bigIntToBytes(mBig);
    byte[] result = UInt256Algo.addMod(a, b, m);
    byte[] expected = bigIntToBytes(aBig.add(bBig).mod(mBig));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testAddMod_resultSmallerThanModulus() {
    byte[] a = new byte[] {3};
    byte[] b = new byte[] {4};
    byte[] m = new byte[] {10};
    byte[] result = UInt256Algo.addMod(a, b, m);
    byte[] expected = new byte[] {7};
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
    // Verify result < modulus
    assertThat(UInt256Algo.compare(result, m)).isLessThan(0);
  }

  // endregion

  // region Modular Multiplication Tests

  @Test
  public void testMulMod_simple() {
    byte[] a = new byte[] {5};
    byte[] b = new byte[] {7};
    byte[] m = new byte[] {11};
    byte[] result = UInt256Algo.mulMod(a, b, m);
    byte[] expected = new byte[] {2}; // (5 * 7) % 11 = 35 % 11 = 2
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testMulMod_withZeroMultiplicand() {
    byte[] zero = new byte[0];
    byte[] b = new byte[] {7};
    byte[] m = new byte[] {11};
    byte[] result = UInt256Algo.mulMod(zero, b, m);
    assertArrayEquals(new byte[0], result);
  }

  @Test
  public void testMulMod_withZeroModulus() {
    byte[] a = new byte[] {5};
    byte[] b = new byte[] {7};
    byte[] zero = new byte[0];
    byte[] result = UInt256Algo.mulMod(a, b, zero);
    assertArrayEquals(new byte[0], result);
  }

  @Test
  public void testMulMod_withOne() {
    byte[] a = new byte[] {42};
    byte[] one = new byte[] {1};
    byte[] m = new byte[] {11};
    byte[] result = UInt256Algo.mulMod(a, one, m);
    byte[] expected = new byte[] {9}; // 42 % 11 = 9
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testMulMod_fromBytes() {
    byte[] a = new byte[] {0, 0, 0, 0, 0, 0, 0, 1};
    byte[] b = new byte[] {1};
    byte[] m = new byte[] {(byte) 0xFF};
    byte[] result = UInt256Algo.mulMod(a, b, m);
    byte[] expected = new byte[] {1};
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testMulMod_largeNumbers() {
    BigInteger aBig = new BigInteger("123456789abcdef", 16);
    BigInteger bBig = new BigInteger("fedcba987654321", 16);
    BigInteger mBig = new BigInteger("1000000007", 16);
    byte[] a = bigIntToBytes(aBig);
    byte[] b = bigIntToBytes(bBig);
    byte[] m = bigIntToBytes(mBig);
    byte[] result = UInt256Algo.mulMod(a, b, m);
    byte[] expected = bigIntToBytes(aBig.multiply(bBig).mod(mBig));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  // endregion

  // region Bitwise AND Tests

  @Test
  public void testAnd_simple() {
    byte[] a = new byte[32];
    byte[] b = new byte[32];
    a[31] = (byte) 0b11110000;
    b[31] = (byte) 0b10101010;
    byte[] result = UInt256Algo.and(a, b);
    byte[] expected = new byte[32];
    expected[31] = (byte) 0b10100000;
    assertArrayEquals(expected, result);
  }

  @Test
  public void testAnd_allZeros() {
    byte[] a = new byte[32];
    byte[] b = new byte[32];
    Arrays.fill(a, (byte) 0xFF);
    byte[] result = UInt256Algo.and(a, b);
    assertArrayEquals(b, result);
  }

  @Test
  public void testAnd_allOnes() {
    byte[] a = new byte[32];
    byte[] b = new byte[32];
    Arrays.fill(a, (byte) 0xFF);
    Arrays.fill(b, (byte) 0xFF);
    byte[] result = UInt256Algo.and(a, b);
    assertArrayEquals(a, result);
  }

  @Test
  public void testAnd_identity() {
    byte[] a = new byte[32];
    Arrays.fill(a, (byte) 0xAA);
    byte[] allOnes = new byte[32];
    Arrays.fill(allOnes, (byte) 0xFF);
    byte[] result = UInt256Algo.and(a, allOnes);
    assertArrayEquals(a, result);
  }

  // endregion

  // region Bitwise XOR Tests

  @Test
  public void testXor_simple() {
    byte[] a = new byte[32];
    byte[] b = new byte[32];
    a[31] = (byte) 0b11110000;
    b[31] = (byte) 0b10101010;
    byte[] result = UInt256Algo.xor(a, b);
    byte[] expected = new byte[32];
    expected[31] = (byte) 0b01011010;
    assertArrayEquals(expected, result);
  }

  @Test
  public void testXor_withZero() {
    byte[] a = new byte[32];
    Arrays.fill(a, (byte) 0xAA);
    byte[] zero = new byte[32];
    byte[] result = UInt256Algo.xor(a, zero);
    assertArrayEquals(a, result);
  }

  @Test
  public void testXor_withSelf() {
    byte[] a = new byte[32];
    Arrays.fill(a, (byte) 0xAA);
    byte[] result = UInt256Algo.xor(a, a);
    byte[] zero = new byte[32];
    assertArrayEquals(zero, result);
  }

  // endregion

  // region Bitwise OR Tests

  @Test
  public void testOr_simple() {
    byte[] a = new byte[32];
    byte[] b = new byte[32];
    a[31] = (byte) 0b11110000;
    b[31] = (byte) 0b10101010;
    byte[] result = UInt256Algo.or(a, b);
    byte[] expected = new byte[32];
    expected[31] = (byte) 0b11111010;
    assertArrayEquals(expected, result);
  }

  @Test
  public void testOr_withZero() {
    byte[] a = new byte[32];
    Arrays.fill(a, (byte) 0xAA);
    byte[] zero = new byte[32];
    byte[] result = UInt256Algo.or(a, zero);
    assertArrayEquals(a, result);
  }

  @Test
  public void testOr_withSelf() {
    byte[] a = new byte[32];
    Arrays.fill(a, (byte) 0xAA);
    byte[] result = UInt256Algo.or(a, a);
    assertArrayEquals(a, result);
  }

  // endregion

  // region Bitwise NOT Tests

  @Test
  public void testNot_simple() {
    byte[] a = new byte[32];
    a[31] = (byte) 0b11110000;
    byte[] result = UInt256Algo.not(a);
    byte[] expected = new byte[32];
    Arrays.fill(expected, (byte) 0xFF);
    expected[31] = (byte) 0b00001111;
    assertArrayEquals(expected, result);
  }

  @Test
  public void testNot_allZeros() {
    byte[] a = new byte[32];
    byte[] result = UInt256Algo.not(a);
    byte[] expected = new byte[32];
    Arrays.fill(expected, (byte) 0xFF);
    assertArrayEquals(expected, result);
  }

  @Test
  public void testNot_allOnes() {
    byte[] a = new byte[32];
    Arrays.fill(a, (byte) 0xFF);
    byte[] result = UInt256Algo.not(a);
    byte[] expected = new byte[32];
    assertArrayEquals(expected, result);
  }

  @Test
  public void testNot_involutive() {
    byte[] a = new byte[32];
    a[0] = (byte) 0x12;
    a[15] = (byte) 0xAB;
    a[31] = (byte) 0xCD;
    byte[] result = UInt256Algo.not(UInt256Algo.not(a));
    assertArrayEquals(a, result);
  }

  // endregion

  // region Comparison Tests

  @Test
  public void testCompare_equal() {
    byte[] a = new byte[] {1, 2, 3};
    byte[] b = new byte[] {1, 2, 3};
    assertThat(UInt256Algo.compare(a, b)).isEqualTo(0);
  }

  @Test
  public void testCompare_lessThan() {
    byte[] a = new byte[] {1, 2, 3};
    byte[] b = new byte[] {1, 2, 4};
    assertThat(UInt256Algo.compare(a, b)).isLessThan(0);
  }

  @Test
  public void testCompare_greaterThan() {
    byte[] a = new byte[] {1, 2, 4};
    byte[] b = new byte[] {1, 2, 3};
    assertThat(UInt256Algo.compare(a, b)).isGreaterThan(0);
  }

  @Test
  public void testCompare_differentLengths() {
    byte[] a = new byte[] {0, 0, 5};
    byte[] b = new byte[] {5};
    assertThat(UInt256Algo.compare(a, b)).isEqualTo(0);
  }

  @Test
  public void testCompare_leadingZeros() {
    byte[] a = new byte[] {0, 0, 0, 10};
    byte[] b = new byte[] {10};
    assertThat(UInt256Algo.compare(a, b)).isEqualTo(0);
  }

  @Test
  public void testCompare_emptyArrays() {
    byte[] a = new byte[0];
    byte[] b = new byte[0];
    assertThat(UInt256Algo.compare(a, b)).isEqualTo(0);
  }

  @Test
  public void testCompare_emptyVsNonEmpty() {
    byte[] a = new byte[0];
    byte[] b = new byte[] {1};
    assertThat(UInt256Algo.compare(a, b)).isLessThan(0);
    assertThat(UInt256Algo.compare(b, a)).isGreaterThan(0);
  }

  // endregion

  // region Edge Cases and Boundary Tests

  @Test
  public void testMaxValue() {
    byte[] max = new byte[32];
    Arrays.fill(max, (byte) 0xFF);

    // Max AND max = max
    assertArrayEquals(max, UInt256Algo.and(max, max));

    // Max OR max = max
    assertArrayEquals(max, UInt256Algo.or(max, max));

    // Max XOR max = 0
    byte[] zero = new byte[32];
    assertArrayEquals(zero, UInt256Algo.xor(max, max));
  }

  @Test
  public void testMinValue() {
    byte[] zero = new byte[0];

    // 0 + 0 = 0
    assertThat(UInt256Algo.compare(UInt256Algo.addIntWidening(zero, zero), zero)).isEqualTo(0);

    // 0 * anything = 0
    byte[] anything = new byte[] {42};
    assertArrayEquals(new byte[0], UInt256Algo.mul(zero, anything));
  }

  @ParameterizedTest
  @ValueSource(ints = {1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31})
  public void testMod_byPrimes(final int prime) {
    byte[] a = new byte[] {100};
    byte[] p = new byte[] {(byte) prime};
    byte[] result = UInt256Algo.mod(a, p);
    byte[] expected = new byte[] {(byte) (100 % prime)};
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testPowerOfTwo() {
    // Test 2^8 = 256
    byte[] a = new byte[] {1, 0};
    byte[] b = new byte[] {1, 0};
    byte[] result = UInt256Algo.mul(a, b);
    // 256 * 256 = 65536 = 0x010000
    byte[] expected = new byte[] {1, 0, 0};
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Test
  public void testAlternatingPatterns() {
    byte[] pattern1 = new byte[32];
    byte[] pattern2 = new byte[32];
    Arrays.fill(pattern1, (byte) 0xAA); // 10101010
    Arrays.fill(pattern2, (byte) 0x55); // 01010101

    // XOR should give all 1s
    byte[] allOnes = new byte[32];
    Arrays.fill(allOnes, (byte) 0xFF);
    assertArrayEquals(allOnes, UInt256Algo.xor(pattern1, pattern2));

    // AND should give all 0s
    byte[] allZeros = new byte[32];
    assertArrayEquals(allZeros, UInt256Algo.and(pattern1, pattern2));
  }

  // endregion
}
