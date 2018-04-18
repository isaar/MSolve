using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.MaterialsTest
{
    public class ElementWithoutMaterial : IFiniteElement
    {
        private readonly double youngModulus, sectionArea;

        public ElementWithoutMaterial(double youngModulus, double sectionArea)
        {
            this.youngModulus = youngModulus;
            this.sectionArea = sectionArea;
        }

        public double[] GetStiffnessMatrix()
        {
            //... K = E/AIlll.ll
            return new double[0];
        }
    }
}
