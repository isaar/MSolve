using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.MaterialsTest
{
    public class AmbroseMaterialProperty : IMaterialProperty
    {
        public AmbroseMaterialProperty(double youngModulus, double poissonRation, double alpha, double ksi)
        {

        }

        public IMaterialState BuildMaterialState(double[] coordinates)
        {
            return new AmbroseMaterialState(this, coordinates);
        }

        class AmbroseMaterialState : IMaterialState
        {
            internal AmbroseMaterialState(AmbroseMaterialProperty property, double[] coordinates)
            {
            }

            public double[] ConstitutiveMatrix
            {
                get
                {
                    throw new NotImplementedException();
                }
            }

            public bool Modified
            {
                get
                {
                    throw new NotImplementedException();
                }
            }

            public double[] Stresses
            {
                get
                {
                    throw new NotImplementedException();
                }
            }

            public void ClearState()
            {
                throw new NotImplementedException();
            }

            public void ClearStresses()
            {
                throw new NotImplementedException();
            }

            public void ResetModified()
            {
                throw new NotImplementedException();
            }

            public void SaveState()
            {
                throw new NotImplementedException();
            }

            public void UpdateMaterial(double[] strains)
            {
                throw new NotImplementedException();
            }
        }
    }
}
