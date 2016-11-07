using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.Solvers.PSM
{
    public class SolverPSM : IIterativeSolver
    {
        protected readonly Model model;
        protected readonly Dictionary<int, IMatrix2D> preconditionersDictionary, kiiDictionary;
        private readonly Dictionary<int, ILinearSystem> subdomainsDictionary;
        private readonly Dictionary<int, List<int>> boundaryDOFsDictionary, internalDOFsDictionary;

        public int CurrentIteration
        {
            get { throw new NotImplementedException(); }
        }

        public SolverPSM(Model model)
        {
            this.model = model;
            preconditionersDictionary = new Dictionary<int, IMatrix2D>(model.SubdomainsDictionary.Count);
            kiiDictionary = new Dictionary<int, IMatrix2D>(model.SubdomainsDictionary.Count);
            subdomainsDictionary = new Dictionary<int, ILinearSystem>(model.SubdomainsDictionary.Count);
            boundaryDOFsDictionary = new Dictionary<int, List<int>>(model.SubdomainsDictionary.Count);
            internalDOFsDictionary = new Dictionary<int, List<int>>(model.SubdomainsDictionary.Count);
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                subdomainsDictionary.Add(subdomain.ID, new SubdomainSkyline(subdomain));
                boundaryDOFsDictionary.Add(subdomain.ID, new List<int>());
            }
        }

        private void EnumerateBoundaryNodes()
        {
            foreach (Node node in model.NodesDictionary.Values)
            {
                int i = 0;
                foreach (int leftSubdomainID in node.SubdomainsDictionary.Keys)
                {
                    int j = -1;
                    foreach (int rightSubdomainID in node.SubdomainsDictionary.Keys)
                    {
                        j++;
                        if (j <= i) continue;

                        Dictionary<DOFType, int> leftDOFs = model.Subdomains[leftSubdomainID].NodalDOFsDictionary[node.ID];
                        Dictionary<DOFType, int> rightDOFs = model.Subdomains[rightSubdomainID].NodalDOFsDictionary[node.ID];
                        foreach (DOFType dofType in leftDOFs.Keys)
                        {
                            int leftDOF = leftDOFs[dofType];
                            int rightDOF = rightDOFs[dofType];
                            if (leftDOF >= 0 && rightDOF >= 0)
                            {
                                if (boundaryDOFsDictionary[leftSubdomainID].IndexOf(leftDOF) < 0)
                                    boundaryDOFsDictionary[leftSubdomainID].Add(leftDOF);
                                if (boundaryDOFsDictionary[rightSubdomainID].IndexOf(rightDOF) < 0)
                                    boundaryDOFsDictionary[rightSubdomainID].Add(rightDOF);
                            }
                        }
                    }
                    i++;
                }
            }

            foreach (var subdomain in subdomainsDictionary)
                internalDOFsDictionary.Add(subdomain.Key, model.Subdomains[subdomain.Key].NodalDOFsDictionary.SelectMany(x => x.Value.Values).Except(boundaryDOFsDictionary[subdomain.Key]).OrderBy(x => x).ToList<int>());
        }

        private void MakeKiiDictionary()
        {
            foreach (var subdomain in subdomainsDictionary)
            {
                SubdomainSkyline s = (SubdomainSkyline)subdomain.Value;
                SkylineMatrix2D k = (SkylineMatrix2D)s.Matrix;
                //var internalDOFs = model.Subdomains[s.ID].NodalDOFsDictionary.SelectMany(x => x.Value.Values).Except(boundaryDOFsDictionary[s.ID]).OrderBy(x => x).ToArray<int>();
                var internalDOFs = internalDOFsDictionary[s.ID];
                int[] kiiIx = new int[internalDOFs.Count + 1];
                int curIx = 0;
                for (int i = 0; i < internalDOFs.Count; i++)
                {
                    kiiIx[i] = curIx;
                    int row = internalDOFs[i];
                    int fromCol = row - k.RowIndex[row + 1] + k.RowIndex[row] + 1;
                    var newCols = internalDOFs.Count(x => x >= fromCol && x <= row);
                    curIx += newCols;
                }
                kiiIx[internalDOFs.Count] = curIx;

                SkylineMatrix2D kii = new SkylineMatrix2D(kiiIx);
                curIx = 0;
                for (int i = 0; i < internalDOFs.Count; i++)
                {
                    int row = internalDOFs[i];
                    int fromCol = row - k.RowIndex[row + 1] + k.RowIndex[row] + 1;
                    var newCols = internalDOFs.Where(x => x >= fromCol && x <= row).OrderByDescending(x => x).ToArray<int>();
                    for (int j = 0; j < newCols.Length; j++)
                    {
                        kii.Data[curIx] = k.Data[k.RowIndex[row] + row - newCols[j]];
                        curIx++;
                    }
                }

                kii.Factorize(1e-8, new List<IVector>(), new List<int>());
                kiiDictionary.Add(s.ID, kii);
            }
        }

        private void MultiplyKib(int subID, double[] vIn, double[] vOut)
        {
            SkylineMatrix2D k = (SkylineMatrix2D)subdomainsDictionary[subID].Matrix;
            var outputDOFs = internalDOFsDictionary[subID];
            Array.Clear(vOut, 0, outputDOFs.Count);
            int pos = 0;
            for (int i = 0; i < k.Rows; i++)
            {
                int height = k.RowIndex[i + 1] - k.RowIndex[i];
                if (height <= 0) continue;
                int iiIx = outputDOFs.IndexOf(i);
                if (iiIx < 0)
                {
                    pos += height;
                    continue;
                }

                vOut[iiIx] += k.Data[pos] * vIn[i];
                pos++;
                for (int j = 0; j < height - 1; j++)
                {
                    int row = i - j - 1;
                    int iiIxRow = outputDOFs.IndexOf(row);
                    if (iiIxRow > -1)
                        vOut[iiIxRow] += k.Data[pos] * vIn[i];
                    vOut[iiIx] += k.Data[pos] * vIn[row];
                    pos++;
                }
            }
        }

        private void MultiplyKbi(int subID, double[] vIn, double[] vOut)
        {
            SkylineMatrix2D k = (SkylineMatrix2D)subdomainsDictionary[subID].Matrix;
            var outputDOFs = boundaryDOFsDictionary[subID];
            Array.Clear(vOut, 0, outputDOFs.Count);
            int pos = 0;
            for (int i = 0; i < k.Rows; i++)
            {
                int height = k.RowIndex[i + 1] - k.RowIndex[i];
                if (height <= 0) continue;
                int iiIx = outputDOFs.IndexOf(i);
                if (iiIx < 0)
                {
                    pos += height;
                    continue;
                }

                vOut[iiIx] += k.Data[pos] * vIn[i];
                pos++;
                for (int j = 0; j < height - 1; j++)
                {
                    int row = i - j - 1;
                    int iiIxRow = outputDOFs.IndexOf(row);
                    if (iiIxRow > -1)
                        vOut[iiIxRow] += k.Data[pos] * vIn[i];
                    vOut[iiIx] += k.Data[pos] * vIn[row];
                    pos++;
                }
            }
        }

        public void Initialize(IVector x, IVector residual, double detf)
        {
            throw new NotImplementedException();
        }

        public void Solve(int maxIterations, double tolerance)
        {
            throw new NotImplementedException();
        }

        public void Initialize()
        {
            throw new NotImplementedException();
        }

        public void Solve()
        {
            throw new NotImplementedException();
        }
    }
}
