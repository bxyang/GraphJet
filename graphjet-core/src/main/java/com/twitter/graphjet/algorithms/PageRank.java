package com.twitter.graphjet.algorithms;

import com.google.common.util.concurrent.AtomicDouble;
import com.twitter.graphjet.bipartite.api.EdgeIterator;
import com.twitter.graphjet.directed.api.OutIndexedDirectedGraph;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.longs.LongIterator;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;

public class PageRank {
  final private OutIndexedDirectedGraph graph;
  final private LongOpenHashSet vertices;
  final private long maxVertexId;
  final private double dampingFactor;
  final private int vertexCount;
  final private double tolerance;

  public PageRank(OutIndexedDirectedGraph graph, LongOpenHashSet vertices,
                  long maxVertexId, double dampingFactor, double tolerance) {
    this.graph = graph;
    this.vertices = vertices;
    this.maxVertexId = maxVertexId;
    this.dampingFactor = dampingFactor;
    this.vertexCount = vertices.size();
    this.tolerance = tolerance;
  }

  public double deltaOfArrays(double a[], double b[]) {
    double ret = 0.0;
    for (int i = 0; i < a.length; ++i) {
      ret += Math.abs(a[i] - b[i]);
    }
    return ret;
  }

  public double[] iterate(double[] prevPR, double dampingAmount, LongArrayList noOuts, AtomicDouble error) {
    double afterPR[] = new double[(int) (maxVertexId + 1)];
    AtomicDouble dangleSum = new AtomicDouble();
    LongIterator iter = noOuts.iterator();
    while (iter.hasNext()) {
      dangleSum.addAndGet(prevPR[(int) iter.nextLong()]);
    }
    dangleSum.set(dangleSum.get() / vertexCount);

    iter = vertices.iterator();
    while (iter.hasNext()) {
      long v = iter.nextLong();
      int outDegree = graph.getOutDegree(v);
      double outWeight = dampingFactor * prevPR[(int) v] / outDegree;
      EdgeIterator edges = graph.getOutEdges(v);
      while (edges.hasNext()) {
        int nbr = (int) edges.nextLong();
        afterPR[nbr] += outWeight;
      }
      afterPR[(int) v] += dampingAmount + dangleSum.get();
    }

    error.set(deltaOfArrays(prevPR, afterPR));
    return afterPR;
  }

  public double[] initializePR(int nodeCount) {
    double prevPR[] = new double[(int) (maxVertexId + 1)];
    vertices.forEach(v -> prevPR[(int) (long) v] = 1.0 / nodeCount);
    return prevPR;
  }

  public double[] run(int maxIteration) {
    LongArrayList noOuts = new LongArrayList();
    LongIterator iter = noOuts.iterator();
    while (iter.hasNext()) {
      long v = iter.nextLong();
      if (graph.getOutDegree(v) == 0) {
        noOuts.add(v);
      }
    }
    double dampingAmount = (1.0 - dampingFactor) / vertexCount;
    double prevPR[] = initializePR(vertexCount);

    int i = 0;
    AtomicDouble error = new AtomicDouble(Double.MAX_VALUE);
    while (i < maxIteration && error.get() > tolerance) {
      prevPR = iterate(prevPR, dampingAmount, noOuts, error);
      i++;
    }

    return prevPR;
  }
}