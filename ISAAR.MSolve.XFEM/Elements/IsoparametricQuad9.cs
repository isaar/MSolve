using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Elements
{
    class IsoparametricQuad9: ContinuumElement2D
    {
        /// <summary>
        /// The caller assumes responsibility for the gauss points and materials
        /// </summary>
        public IsoparametricQuad9(Node2D[] nodes,
            IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materialsOfGaussPoints): 
            base(nodes, IsoparametricInterpolation2D.Quad4, materialsOfGaussPoints)
        {
        }

        public IsoparametricQuad9(Node2D[] nodes, IFiniteElementMaterial2D commonMaterial) :
            base(nodes, IsoparametricInterpolation2D.Quad4, 
                GaussQuadrature2D.Order2x2.GenerateIntegrationPoints(), commonMaterial)
        {
        }
    }
}
