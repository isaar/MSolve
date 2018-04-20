using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IPolynomialChaosCoefficients
    {
        int PsiSize { get; }
        int[][] PsiBasis { get; }
        int[] PsiSquareNorm { get; }
        double[][] HermitePolynomials { get; }
        List<Tuple<int, int, double>>[] GaussianCoefficients { get; }
        Dictionary<Tuple<int, int>, double>[] LognormalCoefficients { get; }
    }
}
