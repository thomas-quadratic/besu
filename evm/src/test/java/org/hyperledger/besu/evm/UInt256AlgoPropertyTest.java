/*
 * Copyright contributors to Besu.
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

import java.math.BigInteger;
import java.util.Arrays;

import net.jqwik.api.Arbitraries;
import net.jqwik.api.Arbitrary;
import net.jqwik.api.ForAll;
import net.jqwik.api.Property;
import net.jqwik.api.Provide;

/** Property-based tests for UInt256Algo public API */
public class UInt256AlgoPropertyTest {

  // region Test Data Providers

  @Provide
  Arbitrary<byte[]> unsigned1to32() {
    return Arbitraries.bytes()
        .array(byte[].class)
        .ofMinSize(1)
        .ofMaxSize(32)
        .map(UInt256AlgoPropertyTest::clampUnsigned32);
  }

  @Provide
  Arbitrary<byte[]> unsigned0to64() {
    return Arbitraries.bytes()
        .array(byte[].class)
        .ofMinSize(0)
        .ofMaxSize(64)
        .map(UInt256AlgoPropertyTest::clampUnsigned32);
  }

  @Provide
  Arbitrary<byte[]> singleLimbUnsigned1to4() {
    return Arbitraries.bytes()
        .array(byte[].class)
        .ofMinSize(1)
        .ofMaxSize(4)
        .map(UInt256AlgoPropertyTest::clampUnsigned32);
  }

  // endregion

  // region Addition Tests (addA, addB, addC, addD)

  @Property
  void property_addA_matchesBigInteger(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    BigInteger A = toBigUnsigned(a);
    BigInteger B = toBigUnsigned(b);

    // Act
    byte[] result = UInt256Algo.addA(a, b);

    // Assert - compare with BigInteger addition (mod 2^256)
    byte[] expected = bigUnsignedToBytes(A.add(B));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Property
  void property_addB_matchesBigInteger(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    BigInteger A = toBigUnsigned(a);
    BigInteger B = toBigUnsigned(b);

    // Act
    byte[] result = UInt256Algo.addB(a, b);

    // Assert
    byte[] expected = bigUnsignedToBytes(A.add(B));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Property
  void property_addC_matchesBigInteger(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    BigInteger A = toBigUnsigned(a);
    BigInteger B = toBigUnsigned(b);

    // Act
    byte[] result = UInt256Algo.addC(a, b);

    // Assert
    byte[] expected = bigUnsignedToBytes(A.add(B));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Property
  void property_addD_matchesBigInteger(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    BigInteger A = toBigUnsigned(a);
    BigInteger B = toBigUnsigned(b);

    // Act
    byte[] result = UInt256Algo.addD(a, b);

    // Assert
    byte[] expected = bigUnsignedToBytes(A.add(B));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Property
  void property_add_implementations_consistent(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Act
    byte[] resultA = UInt256Algo.addA(a, b);
    byte[] resultB = UInt256Algo.addB(a, b);
    byte[] resultC = UInt256Algo.addC(a, b);
    byte[] resultD = UInt256Algo.addD(a, b);

    // Assert - all implementations should produce same result
    assertThat(UInt256Algo.compare(resultA, resultB)).isEqualTo(0);
    assertThat(UInt256Algo.compare(resultA, resultC)).isEqualTo(0);
    assertThat(UInt256Algo.compare(resultA, resultD)).isEqualTo(0);
  }

  @Property
  void property_add_commutative(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Act & Assert - A + B = B + A
    assertThat(UInt256Algo.compare(UInt256Algo.addA(a, b), UInt256Algo.addA(b, a))).isEqualTo(0);
    assertThat(UInt256Algo.compare(UInt256Algo.addB(a, b), UInt256Algo.addB(b, a))).isEqualTo(0);
    assertThat(UInt256Algo.compare(UInt256Algo.addC(a, b), UInt256Algo.addC(b, a))).isEqualTo(0);
    assertThat(UInt256Algo.compare(UInt256Algo.addD(a, b), UInt256Algo.addD(b, a))).isEqualTo(0);
  }

  @Property
  void property_add_associative(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] c) {
    // Act & Assert - (A + B) + C = A + (B + C) with wrapping
    byte[] left = UInt256Algo.addA(UInt256Algo.addA(a, b), c);
    byte[] right = UInt256Algo.addA(a, UInt256Algo.addA(b, c));
    assertThat(UInt256Algo.compare(left, right)).isEqualTo(0);
  }

  @Property
  void property_add_identity(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] zero = new byte[0];

    // Act & Assert - A + 0 = A
    assertThat(UInt256Algo.compare(UInt256Algo.addA(a, zero), a)).isEqualTo(0);
    assertThat(UInt256Algo.compare(UInt256Algo.addA(zero, a), a)).isEqualTo(0);
  }

  @Property
  void property_add_wrapsAround() {
    // Arrange - max value (2^256 - 1) + 1 should wrap to 0
    byte[] maxValue = new byte[32];
    Arrays.fill(maxValue, (byte) 0xFF);
    byte[] one = new byte[] {1};

    // Act
    byte[] result = UInt256Algo.addA(maxValue, one);

    // Assert - wraps to zero (modulo 2^256)
    byte[] expected = bigUnsignedToBytes(BigInteger.ZERO);
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  // endregion

  // region Multiplication Tests

  @Property
  void property_mul_matchesBigInteger(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    BigInteger A = toBigUnsigned(a);
    BigInteger B = toBigUnsigned(b);

    // Act
    byte[] result = UInt256Algo.mul(a, b);

    // Assert - compare with BigInteger multiplication (mod 2^256)
    byte[] expected = bigUnsignedToBytes(A.multiply(B));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Property
  void property_mul_commutative(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Act & Assert - A * B = B * A
    assertThat(UInt256Algo.compare(UInt256Algo.mul(a, b), UInt256Algo.mul(b, a))).isEqualTo(0);
  }

  @Property
  void property_mul_associative(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] c) {
    // Act & Assert - (A * B) * C = A * (B * C) with wrapping
    byte[] left = UInt256Algo.mul(UInt256Algo.mul(a, b), c);
    byte[] right = UInt256Algo.mul(a, UInt256Algo.mul(b, c));
    assertThat(UInt256Algo.compare(left, right)).isEqualTo(0);
  }

  @Property
  void property_mul_identity(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] one = new byte[] {1};

    // Act & Assert - A * 1 = A
    assertThat(UInt256Algo.compare(UInt256Algo.mul(a, one), a)).isEqualTo(0);
    assertThat(UInt256Algo.compare(UInt256Algo.mul(one, a), a)).isEqualTo(0);
  }

  @Property
  void property_mul_zero(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] zero = new byte[0];

    // Act & Assert - A * 0 = 0
    assertThat(UInt256Algo.mul(a, zero)).isEmpty();
    assertThat(UInt256Algo.mul(zero, a)).isEmpty();
  }

  @Property
  void property_mul_distributive(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] c) {
    // Act - A * (B + C) = (A * B) + (A * C) with wrapping
    byte[] left = UInt256Algo.mul(a, UInt256Algo.addA(b, c));
    byte[] right = UInt256Algo.addA(UInt256Algo.mul(a, b), UInt256Algo.mul(a, c));

    // Assert
    assertThat(UInt256Algo.compare(left, right)).isEqualTo(0);
  }

  // endregion

  // region Modulo Tests

  @Property
  void property_mod_matchesBigInteger(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] m) {
    // Arrange
    BigInteger A = toBigUnsigned(a);
    BigInteger M = toBigUnsigned(m);

    // Act
    byte[] result = UInt256Algo.mod(a, m);

    // Assert
    byte[] expected = (M.signum() == 0) ? new byte[0] : bigUnsignedToBytes(A.mod(M));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Property
  void property_mod_singleLimb_matchesBigInteger(
      @ForAll("singleLimbUnsigned1to4") final byte[] a,
      @ForAll("singleLimbUnsigned1to4") final byte[] m) {
    // Arrange
    BigInteger A = toBigUnsigned(a);
    BigInteger M = toBigUnsigned(m);

    // Act
    byte[] result = UInt256Algo.mod(a, m);

    // Assert
    byte[] expected = (M.signum() == 0) ? new byte[0] : bigUnsignedToBytes(A.mod(M));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Property
  void property_mod_resultSmallerThanModulus(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] m) {
    // Arrange
    BigInteger M = toBigUnsigned(m);
    if (M.signum() == 0) return; // Skip if modulus is zero

    // Act
    byte[] result = UInt256Algo.mod(a, m);
    BigInteger R = toBigUnsigned(result);

    // Assert - result < modulus
    assertThat(R.compareTo(M)).isLessThan(0);
  }

  @Property
  void property_mod_byZero(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] zero = new byte[0];

    // Act & Assert - mod by zero returns empty array
    assertThat(UInt256Algo.mod(a, zero)).isEmpty();
  }

  @Property
  void property_mod_zeroModAnything(@ForAll("unsigned1to32") final byte[] m) {
    // Arrange
    byte[] zero = new byte[0];

    // Act
    byte[] result = UInt256Algo.mod(zero, m);

    // Assert - 0 mod M = 0 (or empty if M is zero)
    BigInteger M = toBigUnsigned(m);
    if (M.signum() == 0) {
      assertThat(result).isEmpty();
    } else {
      assertThat(result).isEmpty();
    }
  }

  @Property
  void property_mod_idempotent(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] m) {
    // Act - (A mod M) mod M = A mod M
    byte[] once = UInt256Algo.mod(a, m);
    byte[] twice = UInt256Algo.mod(once, m);

    // Assert
    assertThat(UInt256Algo.compare(once, twice)).isEqualTo(0);
  }

  // endregion

  // region Signed Modulo Tests

  @Property
  void property_signedMod_matchesEvmSemantics(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] m) {
    // Arrange - pad to 32 bytes for signed interpretation
    byte[] a32 = padTo32Bytes(a);
    byte[] m32 = padTo32Bytes(m);
    BigInteger A = new BigInteger(a32);
    BigInteger M = new BigInteger(m32);

    // Act
    byte[] result = UInt256Algo.signedMod(a32, m32);

    // Assert
    byte[] expected = (M.signum() == 0) ? new byte[0] : computeSignedModExpected(A, M);
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Property
  void property_signedMod_byZero(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] zero = new byte[0];

    // Act & Assert
    assertThat(UInt256Algo.signedMod(a, zero)).isEmpty();
  }

  // endregion

  // region Modular Addition Tests

  @Property
  void property_addMod_matchesBigInteger(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] m) {
    // Arrange
    BigInteger A = toBigUnsigned(a);
    BigInteger B = toBigUnsigned(b);
    BigInteger M = toBigUnsigned(m);

    // Act
    byte[] result = UInt256Algo.addMod(a, b, m);

    // Assert
    byte[] expected = (M.signum() == 0) ? new byte[0] : bigUnsignedToBytes(A.add(B).mod(M));
    System.out.println(String.format("%s + %s (mod %s)", Arrays.toString(a), Arrays.toString(b), Arrays.toString(m)));
    System.out.println(String.format("Expected %s", Arrays.toString(expected)));
    System.out.println(String.format("Was %s", Arrays.toString(result)));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Property
  void property_addMod_commutative(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] m) {
    // Act & Assert - (A + B) mod M = (B + A) mod M
    assertThat(UInt256Algo.compare(UInt256Algo.addMod(a, b, m), UInt256Algo.addMod(b, a, m)))
        .isEqualTo(0);
  }

  @Property
  void property_addMod_associative(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] c,
      @ForAll("unsigned1to32") final byte[] m) {
    // Act & Assert - ((A + B) + C) mod M = (A + (B + C)) mod M
    byte[] left = UInt256Algo.addMod(UInt256Algo.addMod(a, b, m), c, m);
    byte[] right = UInt256Algo.addMod(a, UInt256Algo.addMod(b, c, m), m);
    assertThat(UInt256Algo.compare(left, right)).isEqualTo(0);
  }

  @Property
  void property_addMod_identity(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] m) {
    // Arrange
    byte[] zero = new byte[0];

    // Act & Assert - (A + 0) mod M = A mod M
    assertThat(UInt256Algo.compare(UInt256Algo.addMod(a, zero, m), UInt256Algo.mod(a, m)))
        .isEqualTo(0);
    assertThat(UInt256Algo.compare(UInt256Algo.addMod(zero, a, m), UInt256Algo.mod(a, m)))
        .isEqualTo(0);
  }

  @Property
  void property_addMod_byZeroModulus(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] zero = new byte[0];

    // Act & Assert - (A + B) mod 0 = empty
    assertThat(UInt256Algo.addMod(a, b, zero)).isEmpty();
  }

  @Property
  void property_addMod_resultSmallerThanModulus(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] m) {
    // Arrange
    BigInteger M = toBigUnsigned(m);
    if (M.signum() == 0 || M.equals(BigInteger.ONE)) return; // Skip if modulus is 0 or 1

    // Act
    byte[] result = UInt256Algo.addMod(a, b, m);
    BigInteger R = toBigUnsigned(result);

    // Assert - result < modulus
    assertThat(R.compareTo(M)).isLessThan(0);
  }

  // endregion

  // region Modular Multiplication Tests

  @Property
  void property_mulMod_matchesBigInteger(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] m) {
    // Arrange
    BigInteger A = toBigUnsigned(a);
    BigInteger B = toBigUnsigned(b);
    BigInteger M = toBigUnsigned(m);

    // Act
    byte[] result = UInt256Algo.mulMod(a, b, m);

    // Assert
    byte[] expected = (M.signum() == 0) ? new byte[0] : bigUnsignedToBytes(A.multiply(B).mod(M));
    assertThat(UInt256Algo.compare(result, expected)).isEqualTo(0);
  }

  @Property
  void property_mulMod_commutative(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] m) {
    // Act & Assert - (A * B) mod M = (B * A) mod M
    assertThat(UInt256Algo.compare(UInt256Algo.mulMod(a, b, m), UInt256Algo.mulMod(b, a, m)))
        .isEqualTo(0);
  }

  @Property
  void property_mulMod_associative(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] c,
      @ForAll("unsigned1to32") final byte[] m) {
    // Act & Assert - ((A * B) * C) mod M = (A * (B * C)) mod M
    byte[] left = UInt256Algo.mulMod(UInt256Algo.mulMod(a, b, m), c, m);
    byte[] right = UInt256Algo.mulMod(a, UInt256Algo.mulMod(b, c, m), m);
    assertThat(UInt256Algo.compare(left, right)).isEqualTo(0);
  }

  @Property
  void property_mulMod_identity(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] m) {
    // Arrange
    byte[] one = new byte[] {1};

    // Act & Assert - (A * 1) mod M = A mod M
    assertThat(UInt256Algo.compare(UInt256Algo.mulMod(a, one, m), UInt256Algo.mod(a, m)))
        .isEqualTo(0);
    assertThat(UInt256Algo.compare(UInt256Algo.mulMod(one, a, m), UInt256Algo.mod(a, m)))
        .isEqualTo(0);
  }

  @Property
  void property_mulMod_zero(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] m) {
    // Arrange
    byte[] zero = new byte[0];

    // Act & Assert - (A * 0) mod M = 0 (or empty if M is 0)
    assertThat(UInt256Algo.mulMod(a, zero, m)).isEmpty();
    assertThat(UInt256Algo.mulMod(zero, a, m)).isEmpty();
  }

  @Property
  void property_mulMod_byZeroModulus(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] zero = new byte[0];

    // Act & Assert - (A * B) mod 0 = empty
    assertThat(UInt256Algo.mulMod(a, b, zero)).isEmpty();
  }

  @Property
  void property_mulMod_resultSmallerThanModulus(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] m) {
    // Arrange
    BigInteger M = toBigUnsigned(m);
    if (M.signum() == 0 || M.equals(BigInteger.ONE)) return; // Skip if modulus is 0 or 1

    // Act
    byte[] result = UInt256Algo.mulMod(a, b, m);
    BigInteger R = toBigUnsigned(result);

    // Assert - result < modulus
    assertThat(R.compareTo(M)).isLessThan(0);
  }

  @Property
  void property_mulMod_distributive(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] c,
      @ForAll("unsigned1to32") final byte[] m) {
    // Act - (A * (B + C)) mod M = ((A * B) + (A * C)) mod M
    byte[] left = UInt256Algo.mulMod(a, UInt256Algo.addMod(b, c, m), m);
    byte[] right = UInt256Algo.addMod(UInt256Algo.mulMod(a, b, m), UInt256Algo.mulMod(a, c, m), m);

    // Assert
    assertThat(UInt256Algo.compare(left, right)).isEqualTo(0);
  }

  // endregion

  // region Bitwise AND Tests

  @Property
  void property_and_matchesBigInteger(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);
    BigInteger A = toBigUnsigned(a32);
    BigInteger B = toBigUnsigned(b32);

    // Act
    byte[] result = UInt256Algo.and(a32, b32);

    // Assert
    byte[] expected = bigUnsignedToBytes32(A.and(B));
    assertThat(result).containsExactly(expected);
  }

  @Property
  void property_and_commutative(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);

    // Act & Assert - A & B = B & A
    assertThat(UInt256Algo.and(a32, b32)).containsExactly(UInt256Algo.and(b32, a32));
  }

  @Property
  void property_and_associative(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] c) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);
    byte[] c32 = padTo32Bytes(c);

    // Act & Assert - (A & B) & C = A & (B & C)
    assertThat(UInt256Algo.and(UInt256Algo.and(a32, b32), c32))
        .containsExactly(UInt256Algo.and(a32, UInt256Algo.and(b32, c32)));
  }

  @Property
  void property_and_identity(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] allOnes = new byte[32];
    Arrays.fill(allOnes, (byte) 0xFF);

    // Act & Assert - A & 0xFF...FF = A
    assertThat(UInt256Algo.and(a32, allOnes)).containsExactly(a32);
  }

  @Property
  void property_and_zero(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] zero = new byte[32];

    // Act & Assert - A & 0 = 0
    assertThat(UInt256Algo.and(a32, zero)).containsExactly(zero);
  }

  @Property
  void property_and_self(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);

    // Act & Assert - A & A = A (idempotent)
    assertThat(UInt256Algo.and(a32, a32)).containsExactly(a32);
  }

  // endregion

  // region Bitwise XOR Tests

  @Property
  void property_xor_matchesBigInteger(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);
    BigInteger A = toBigUnsigned(a32);
    BigInteger B = toBigUnsigned(b32);

    // Act
    byte[] result = UInt256Algo.xor(a32, b32);

    // Assert
    byte[] expected = bigUnsignedToBytes32(A.xor(B));
    assertThat(result).containsExactly(expected);
  }

  @Property
  void property_xor_commutative(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);

    // Act & Assert - A ^ B = B ^ A
    assertThat(UInt256Algo.xor(a32, b32)).containsExactly(UInt256Algo.xor(b32, a32));
  }

  @Property
  void property_xor_associative(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] c) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);
    byte[] c32 = padTo32Bytes(c);

    // Act & Assert - (A ^ B) ^ C = A ^ (B ^ C)
    assertThat(UInt256Algo.xor(UInt256Algo.xor(a32, b32), c32))
        .containsExactly(UInt256Algo.xor(a32, UInt256Algo.xor(b32, c32)));
  }

  @Property
  void property_xor_identity(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] zero = new byte[32];

    // Act & Assert - A ^ 0 = A
    assertThat(UInt256Algo.xor(a32, zero)).containsExactly(a32);
  }

  @Property
  void property_xor_self(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] zero = new byte[32];

    // Act & Assert - A ^ A = 0 (self-inverse)
    assertThat(UInt256Algo.xor(a32, a32)).containsExactly(zero);
  }

  @Property
  void property_xor_involutive(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);

    // Act & Assert - (A ^ B) ^ B = A
    assertThat(UInt256Algo.xor(UInt256Algo.xor(a32, b32), b32)).containsExactly(a32);
  }

  // endregion

  // region Bitwise OR Tests

  @Property
  void property_or_matchesBigInteger(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);
    BigInteger A = toBigUnsigned(a32);
    BigInteger B = toBigUnsigned(b32);

    // Act
    byte[] result = UInt256Algo.or(a32, b32);

    // Assert
    byte[] expected = bigUnsignedToBytes32(A.or(B));
    assertThat(result).containsExactly(expected);
  }

  @Property
  void property_or_commutative(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);

    // Act & Assert - A | B = B | A
    assertThat(UInt256Algo.or(a32, b32)).containsExactly(UInt256Algo.or(b32, a32));
  }

  @Property
  void property_or_associative(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] c) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);
    byte[] c32 = padTo32Bytes(c);

    // Act & Assert - (A | B) | C = A | (B | C)
    assertThat(UInt256Algo.or(UInt256Algo.or(a32, b32), c32))
        .containsExactly(UInt256Algo.or(a32, UInt256Algo.or(b32, c32)));
  }

  @Property
  void property_or_identity(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] zero = new byte[32];

    // Act & Assert - A | 0 = A
    assertThat(UInt256Algo.or(a32, zero)).containsExactly(a32);
  }

  @Property
  void property_or_self(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);

    // Act & Assert - A | A = A (idempotent)
    assertThat(UInt256Algo.or(a32, a32)).containsExactly(a32);
  }

  @Property
  void property_or_absorbs_and(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);

    // Act & Assert - A | (A & B) = A (absorption law)
    assertThat(UInt256Algo.or(a32, UInt256Algo.and(a32, b32))).containsExactly(a32);
  }

  // endregion

  // region Bitwise NOT Tests

  @Property
  void property_not_matchesBigInteger(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    BigInteger A = toBigUnsigned(a32);
    BigInteger allOnes = BigInteger.ONE.shiftLeft(256).subtract(BigInteger.ONE);

    // Act
    byte[] result = UInt256Algo.not(a32);

    // Assert - ~A = A XOR 0xFF...FF
    byte[] expected = bigUnsignedToBytes32(A.xor(allOnes));
    assertThat(result).containsExactly(expected);
  }

  @Property
  void property_not_involutive(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);

    // Act & Assert - ~~A = A
    assertThat(UInt256Algo.not(UInt256Algo.not(a32))).containsExactly(a32);
  }

  @Property
  void property_not_with_and_is_zero(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] zero = new byte[32];

    // Act & Assert - A & ~A = 0
    assertThat(UInt256Algo.and(a32, UInt256Algo.not(a32))).containsExactly(zero);
  }

  @Property
  void property_not_with_or_is_allOnes(@ForAll("unsigned1to32") final byte[] a) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] allOnes = new byte[32];
    Arrays.fill(allOnes, (byte) 0xFF);

    // Act & Assert - A | ~A = 0xFF...FF
    assertThat(UInt256Algo.or(a32, UInt256Algo.not(a32))).containsExactly(allOnes);
  }

  @Property
  void property_not_de_morgans_and(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);

    // Act & Assert - ~(A & B) = ~A | ~B
    byte[] left = UInt256Algo.not(UInt256Algo.and(a32, b32));
    byte[] right = UInt256Algo.or(UInt256Algo.not(a32), UInt256Algo.not(b32));
    assertThat(left).containsExactly(right);
  }

  @Property
  void property_not_de_morgans_or(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);

    // Act & Assert - ~(A | B) = ~A & ~B
    byte[] left = UInt256Algo.not(UInt256Algo.or(a32, b32));
    byte[] right = UInt256Algo.and(UInt256Algo.not(a32), UInt256Algo.not(b32));
    assertThat(left).containsExactly(right);
  }

  // endregion

  // region Combined Bitwise Properties

  @Property
  void property_and_distributes_over_xor(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] c) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);
    byte[] c32 = padTo32Bytes(c);

    // Act & Assert - A & (B ^ C) = (A & B) ^ (A & C)
    byte[] left = UInt256Algo.and(a32, UInt256Algo.xor(b32, c32));
    byte[] right = UInt256Algo.xor(UInt256Algo.and(a32, b32), UInt256Algo.and(a32, c32));
    assertThat(left).containsExactly(right);
  }

  @Property
  void property_or_distributes_over_and(
      @ForAll("unsigned1to32") final byte[] a,
      @ForAll("unsigned1to32") final byte[] b,
      @ForAll("unsigned1to32") final byte[] c) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);
    byte[] c32 = padTo32Bytes(c);

    // Act & Assert - A | (B & C) = (A | B) & (A | C)
    byte[] left = UInt256Algo.or(a32, UInt256Algo.and(b32, c32));
    byte[] right = UInt256Algo.and(UInt256Algo.or(a32, b32), UInt256Algo.or(a32, c32));
    assertThat(left).containsExactly(right);
  }

  @Property
  void property_or_equals_and_plus_xor(
      @ForAll("unsigned1to32") final byte[] a, @ForAll("unsigned1to32") final byte[] b) {
    // Arrange
    byte[] a32 = padTo32Bytes(a);
    byte[] b32 = padTo32Bytes(b);

    // Act & Assert - A | B = (A & B) | (A ^ B)
    byte[] left = UInt256Algo.or(a32, b32);
    byte[] right = UInt256Algo.or(UInt256Algo.and(a32, b32), UInt256Algo.xor(a32, b32));
    assertThat(left).containsExactly(right);
  }

  // endregion

  // region Utility Methods

  private static byte[] clampUnsigned32(final byte[] any) {
    if (any.length == 0) {
      return new byte[] {0};
    }
    int len = Math.max(1, Math.min(32, any.length));
    byte[] out = new byte[len];
    System.arraycopy(any, 0, out, 0, len);
    return out;
  }

  private static byte[] padTo32Bytes(final byte[] a) {
    if (a.length >= 32) {
      byte[] result = new byte[32];
      System.arraycopy(a, Math.max(0, a.length - 32), result, 0, 32);
      return result;
    }
    byte[] result = new byte[32];
    System.arraycopy(a, 0, result, 32 - a.length, a.length);
    return result;
  }

  private static byte[] bigUnsignedToBytes(final BigInteger x) {
    BigInteger y = x.mod(BigInteger.ONE.shiftLeft(256));

    byte[] ba = y.toByteArray();
    if (ba.length == 0) {
      return new byte[0];
    }

    // Remove sign bit if present
    if (ba.length > 32 && ba[0] == 0) {
      byte[] out = new byte[32];
      System.arraycopy(ba, ba.length - 32, out, 0, 32);
      return out;
    }

    if (ba.length <= 32) {
      // Remove leading zeros
      int firstNonZero = 0;
      while (firstNonZero < ba.length && ba[firstNonZero] == 0) {
        firstNonZero++;
      }
      if (firstNonZero == ba.length) {
        return new byte[0];
      }
      byte[] out = new byte[ba.length - firstNonZero];
      System.arraycopy(ba, firstNonZero, out, 0, out.length);
      return out;
    }

    // If bigger than 32, take lower 32 bytes
    byte[] out = new byte[32];
    System.arraycopy(ba, ba.length - 32, out, 0, 32);
    return out;
  }

  private static byte[] bigUnsignedToBytes32(final BigInteger x) {
    BigInteger y = x.mod(BigInteger.ONE.shiftLeft(256));

    byte[] ba = y.toByteArray();
    if (ba.length == 0) {
      return new byte[32];
    }

    if (ba.length == 32) {
      return ba;
    }

    if (ba.length < 32) {
      byte[] out = new byte[32];
      System.arraycopy(ba, 0, out, 32 - ba.length, ba.length);
      return out;
    }

    // If bigger than 32, take lower 32 bytes
    byte[] out = new byte[32];
    System.arraycopy(ba, ba.length - 32, out, 0, 32);
    return out;
  }

  private static BigInteger toBigUnsigned(final byte[] be) {
    if (be.length == 0) return BigInteger.ZERO;
    return new BigInteger(1, be);
  }

  private static byte[] computeSignedModExpected(final BigInteger A, final BigInteger M) {
    BigInteger r = A.abs().mod(M.abs());

    if (A.signum() < 0 && r.signum() != 0) {
      return padNegative(r);
    }

    return bigUnsignedToBytes(r);
  }

  private static byte[] padNegative(final BigInteger r) {
    BigInteger neg = r.negate();
    byte[] rb = neg.toByteArray();
    if (rb.length >= 32) {
      byte[] result = new byte[32];
      System.arraycopy(rb, Math.max(0, rb.length - 32), result, 0, 32);
      return result;
    }
    byte[] padded = new byte[32];
    Arrays.fill(padded, (byte) 0xFF);
    System.arraycopy(rb, 0, padded, 32 - rb.length, rb.length);
    return padded;
  }

  // endregion
}
