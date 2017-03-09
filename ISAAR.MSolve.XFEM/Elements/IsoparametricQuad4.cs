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
    // I could very easily make ContinuumElement2D generic on the Interpolation type. However, these concrete classes
    // are more readable (and expected in FEM software) and may be useful for geometric operations later on.
    class IsoparametricQuad4: ContinuumElement2D
    {
        /// <summary>
        /// The caller assumes responsibility for the the nodes, gauss points and materials
        /// </summary>
        public IsoparametricQuad4(IReadOnlyList<Node2D> nodes, 
            IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materialsOfGaussPoints): 
            base(nodes, IsoparametricInterpolation2D.Quad4, materialsOfGaussPoints)
        {
        }

        /// <summary>
        /// The caller assumes responsibility for the the nodes.
        /// </summary>
        public IsoparametricQuad4(IReadOnlyList<Node2D> nodes, IFiniteElementMaterial2D commonMaterial) :
            base(nodes, IsoparametricInterpolation2D.Quad4, 
                GaussLegendre2D .Order2x2.GenerateIntegrationPoints(), commonMaterial)
        {
        }

        public IsoparametricQuad4(IReadOnlyList<Node2D> nodes) :
            base(nodes, IsoparametricInterpolation2D.Quad4, null)
        {
        }
    }
}
