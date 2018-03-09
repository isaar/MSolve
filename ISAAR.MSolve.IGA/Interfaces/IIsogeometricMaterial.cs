using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.IGA.Interfaces
{
    public interface IIsogeometricMaterial : ICloneable
    {
        int ID { get; }
        bool Modified { get; }
        void ResetModified();
        double YoungModulus { get; set; }
        double PoissonRatio { get; set; }
        double[] Coordinates { get; set; }
    }
}
