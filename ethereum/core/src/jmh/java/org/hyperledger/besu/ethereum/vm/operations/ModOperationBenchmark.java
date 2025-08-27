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
package org.hyperledger.besu.ethereum.vm.operations;

import org.hyperledger.besu.evm.frame.MessageFrame;
import org.hyperledger.besu.evm.operation.ModOperation;
import org.hyperledger.besu.evm.operation.Operation;

import java.math.BigInteger;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.tuweni.bytes.Bytes;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.infra.Blackhole;

public class ModOperationBenchmark extends BinaryOperationBenchmark {
  // Benches for a % b

  protected Bytes[] aSmallerPool;
  protected Bytes[] bSmallerPool;
  protected Bytes[] aIntPool;
  protected Bytes[] bIntPool;
  protected Bytes[] aLongPool;
  protected Bytes[] bLongPool;
  protected Bytes[] aSmallPool;
  protected Bytes[] bSmallPool;
  protected Bytes[] aBigPool;
  protected Bytes[] bBigPool;

  @Setup
  @Override
  public void setUp() {
    frame = BenchmarkHelper.createMessageFrame();

    // executeOperation Pool
    aPool = new Bytes[SAMPLE_SIZE];
    bPool = new Bytes[SAMPLE_SIZE];

    // Pool for a < b
    aSmallerPool = new Bytes[SAMPLE_SIZE];
    bSmallerPool = new Bytes[SAMPLE_SIZE];

    // Pool for b fits in int
    aIntPool = new Bytes[SAMPLE_SIZE];
    bIntPool = new Bytes[SAMPLE_SIZE];

    // Pool for b fits in long
    aLongPool = new Bytes[SAMPLE_SIZE];
    bLongPool = new Bytes[SAMPLE_SIZE];

    // Pool for small b, but bigger than long
    aSmallPool = new Bytes[SAMPLE_SIZE];
    bSmallPool = new Bytes[SAMPLE_SIZE];

    // Pool for big b, more than half the bytes
    aBigPool = new Bytes[SAMPLE_SIZE];
    bBigPool = new Bytes[SAMPLE_SIZE];

    final ThreadLocalRandom random = ThreadLocalRandom.current();

    for (int i = 0; i < SAMPLE_SIZE; i++) {
      final int bSize = 1 + random.nextInt(32);
      final int aSize = 1 + random.nextInt(32);
      final int bIntSize = 1 + random.nextInt(4);
      final int aIntSize = 5 + random.nextInt(28);
      final int bLongSize = 9 + random.nextInt(8);
      final int aLongSize = 9 + random.nextInt(24);
      final int bSmallSize = 9 + random.nextInt(8);
      final int aSmallSize = bSmallSize + random.nextInt(16);
      final int bBigSize = 17 + random.nextInt(16);
      final int aBigSize = Math.min(32, bBigSize + random.nextInt(16));

      final byte[] a = new byte[aSize];
      final byte[] b = new byte[bSize];
      final byte[] aInt = new byte[aIntSize];
      final byte[] bInt = new byte[bIntSize];
      final byte[] aLong = new byte[aLongSize];
      final byte[] bLong = new byte[bLongSize];
      final byte[] aSmall = new byte[aSmallSize];
      final byte[] bSmall = new byte[bSmallSize];
      final byte[] aBig = new byte[aBigSize];
      final byte[] bBig = new byte[bBigSize];

      random.nextBytes(a);
      random.nextBytes(b);
      random.nextBytes(aInt);
      random.nextBytes(bInt);
      random.nextBytes(aLong);
      random.nextBytes(bLong);
      random.nextBytes(aSmall);
      random.nextBytes(bSmall);
      random.nextBytes(aBig);
      random.nextBytes(bBig);

      BigInteger aNum = new BigInteger(1, a);
      BigInteger bNum = new BigInteger(1, b);

      if (aNum.compareTo(bNum) < 0) {
        aPool[i] = Bytes.wrap(b);
        bPool[i] = Bytes.wrap(a);
        aSmallerPool[i] = Bytes.wrap(a);
        bSmallerPool[i] = Bytes.wrap(b);
      } else {
        aPool[i] = Bytes.wrap(a);
        bPool[i] = Bytes.wrap(b);
        aSmallerPool[i] = Bytes.wrap(b);
        bSmallerPool[i] = Bytes.wrap(a);
      }

      aIntPool[i] = Bytes.wrap(aInt);
      bIntPool[i] = Bytes.wrap(bInt);
      aLongPool[i] = Bytes.wrap(aLong);
      bLongPool[i] = Bytes.wrap(bLong);
      aSmallPool[i] = Bytes.wrap(aSmall);
      bSmallPool[i] = Bytes.wrap(bSmall);
      aBigPool[i] = Bytes.wrap(aBig);
      bBigPool[i] = Bytes.wrap(bBig);
    }

    index = 0;
  }

  @Benchmark
  public void onZeroModulus(final Blackhole blackhole) {
    final int i = index;
    index = (index + 1) % SAMPLE_SIZE;

    frame.pushStackItem(Bytes.EMPTY);
    frame.pushStackItem(aPool[i]);

    blackhole.consume(invoke(frame));

    frame.popStackItem();
  }

  @Benchmark
  public void onLargerModulus(final Blackhole blackhole) {
    final int i = index;
    index = (index + 1) % SAMPLE_SIZE;

    frame.pushStackItem(bSmallerPool[i]);
    frame.pushStackItem(aSmallPool[i]);

    blackhole.consume(invoke(frame));

    frame.popStackItem();
  }

  @Benchmark
  public void onIntModulus(final Blackhole blackhole) {
    final int i = index;
    index = (index + 1) % SAMPLE_SIZE;

    frame.pushStackItem(bIntPool[i]);
    frame.pushStackItem(aIntPool[i]);

    blackhole.consume(invoke(frame));

    frame.popStackItem();
  }

  @Benchmark
  public void onLongModulus(final Blackhole blackhole) {
    final int i = index;
    index = (index + 1) % SAMPLE_SIZE;

    frame.pushStackItem(bLongPool[i]);
    frame.pushStackItem(aLongPool[i]);

    blackhole.consume(invoke(frame));

    frame.popStackItem();
  }

  @Benchmark
  public void onSmallModulus(final Blackhole blackhole) {
    final int i = index;
    index = (index + 1) % SAMPLE_SIZE;

    frame.pushStackItem(bSmallPool[i]);
    frame.pushStackItem(aSmallPool[i]);

    blackhole.consume(invoke(frame));

    frame.popStackItem();
  }

  @Benchmark
  public void onBigModulus(final Blackhole blackhole) {
    final int i = index;
    index = (index + 1) % SAMPLE_SIZE;

    frame.pushStackItem(bBigPool[i]);
    frame.pushStackItem(aBigPool[i]);

    blackhole.consume(invoke(frame));

    frame.popStackItem();
  }

  @Override
  protected Operation.OperationResult invoke(final MessageFrame frame) {
    return ModOperation.staticOperation(frame);
  }
}
