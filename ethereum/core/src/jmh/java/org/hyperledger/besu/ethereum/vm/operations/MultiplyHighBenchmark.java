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

import java.util.concurrent.TimeUnit;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.tuweni.bytes.Bytes;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

@State(Scope.Thread)
@Warmup(iterations = 2, time = 1, timeUnit = TimeUnit.SECONDS)
@OutputTimeUnit(value = TimeUnit.NANOSECONDS)
@Measurement(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@BenchmarkMode(Mode.AverageTime)
public class MultiplyHighBenchmark {
  protected static final int SAMPLE_SIZE = 80_000;
  protected static final int N_OPS = 100;

  protected long[][] pool;
  protected int index;

  @Setup()
  public void setUp() {
    pool = new long[SAMPLE_SIZE][];
    final ThreadLocalRandom random = ThreadLocalRandom.current();

    for (int i = 0; i < SAMPLE_SIZE; i++) {
      pool[i] = new long[2 * N_OPS];
      for (int j = 0; j < 2 * N_OPS; j += 2) {
        pool[i][j] = random.nextLong();
        pool[i][j + 1] = random.nextLong();
      }
    }
    index = 0;
  }

  @Benchmark
  public void baseline(final Blackhole blackhole) {
    long acc = 0;
    for (int i=0; i < 2 * N_OPS; i += 2) {
      acc = acc | pool[index][i] | pool[index][i + 1];
    }
    blackhole.consume(acc);
    index = (index + 1) % SAMPLE_SIZE;
  }

  @Benchmark
  public void defaultImplem(final Blackhole blackhole) {
    long acc = 0;
    long x;
    long y;
    for (int i=0; i < 2 * N_OPS; i += 2) {
      x = pool[index][i];
      y = pool[index][i + 1];
      long x1 = x >> 32;
      long x2 = x & 0xFFFFFFFFL;
      long y1 = y >> 32;
      long y2 = y & 0xFFFFFFFFL;

      long z2 = x2 * y2;
      long t = x1 * y2 + (z2 >>> 32);
      long z1 = t & 0xFFFFFFFFL;
      long z0 = t >> 32;
      z1 += x2 * y1;

      acc |= (x1 * y1 + z0 + (z1 >> 32));
    }
    blackhole.consume(acc);
    index = (index + 1) % SAMPLE_SIZE;
  }

  @Benchmark
  public void executeOperation(final Blackhole blackhole) {
    long acc = 0;
    for (int i=0; i < 2 * N_OPS; i += 2) {
      acc = acc | Math.multiplyHigh(pool[index][i], pool[index][i + 1]);
    }
    blackhole.consume(acc);
    index = (index + 1) % SAMPLE_SIZE;
  }
}
