using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IFiniteElementMaterial : ICloneable
    {
        int ID { get; }
        bool Modified { get; }
        void ResetModified();
        double YoungModulus { get; set; }
        double PoissonRatio { get; set; }
        double[] Coordinates { get; set; }
    }
}
