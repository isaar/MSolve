using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Functions
{
    class SignFunctionOpposite2D: IHeavisideFunction2D
    {
        public SignFunctionOpposite2D()
        {
        }

        public double EvaluateAt(double signedDistance)
        {
            if (signedDistance > 0.0) return -1.0;
            else if (signedDistance < 0.0) return 1.0;
            else return 0.0;
        }

        public EvaluatedFunction2D EvaluateAllAt(double signedDistance)
        {
            var derivatives = Vector2.Create(0.0, 0.0 );
            if (signedDistance > 0) return new EvaluatedFunction2D(-1.0, derivatives);
            else if (signedDistance < 0) return new EvaluatedFunction2D(1.0, derivatives);
            else return new EvaluatedFunction2D(0.0, derivatives);
        }

        public override string ToString()
        {
            return "Heaviside";
        }
    }
}
