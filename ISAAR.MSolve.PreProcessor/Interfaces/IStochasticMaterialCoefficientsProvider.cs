using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.PreProcessor.Interfaces
{
    public interface IStochasticMaterialCoefficientsProvider
    {
        double[] RandomVariables { get; set; }
        double GetCoefficient(double meanValue, double[] coordinates);
    }
}
