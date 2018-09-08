using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Jintegral
{
    /// <summary>
    /// All quantities are with respect to the local cartesian system of the crack tip.
    /// </summary>
    class AuxiliaryStatesTensors
    {
        /// <summary>
        /// The derivatives of the displacement field of an imaginary pure Mode I (opening crack extension) state 
        /// represented in the local cartesian system of the crack tip. The differentation is also w.r.t. the tip local
        /// cartesian coordinates. Matrix dimensions: 2x2.
        /// </summary>
        public Matrix2by2 DisplacementGradientMode1 { get; }

        /// <summary>
        /// The derivatives of the displacement field of an imaginary pure Mode II (sliding crack extension) state 
        /// represented in the local cartesian system of the crack tip. The differentation is also w.r.t. the tip local 
        /// cartesian coordinates. Matrix dimensions: 2x2.
        /// </summary>
        public Matrix2by2 DisplacementGradientMode2 { get; }

        /// <summary>
        /// The strain tensor of an imaginary pure Mode I (opening crack extension) state represented in the local
        /// cartesian system of the crack tip. Tensor dimensions: 2x2.
        /// </summary>
        public Tensor2D StrainTensorMode1 { get; }

        /// <summary>
        /// The strain tensor of an imaginary pure Mode IΙ (sliding crack extension) state represented in the local
        /// cartesian system of the crack tip. Tensor dimensions: 2x2.
        /// </summary>
        public Tensor2D StrainTensorMode2 { get; }

        /// <summary>
        /// The stress tensor of an imaginary pure Mode I (opening crack extension) state represented in the local
        /// cartesian system of the crack tip. Tensor dimensions: 2x2.
        /// </summary>
        public Tensor2D StressTensorMode1 { get; }

        /// <summary>
        /// The stress tensor of an imaginary pure Mode IΙ (sliding crack extension) state represented in the local
        /// cartesian system of the crack tip. Tensor dimensions: 2x2.
        /// </summary>
        public Tensor2D StressTensorMode2 { get; }

        public AuxiliaryStatesTensors(Matrix2by2 displacementGradientMode1, Matrix2by2 displacementGradientMode2,
            Tensor2D strainTensorMode1, Tensor2D strainTensorMode2, 
            Tensor2D stressTensorMode1, Tensor2D stressTensorMode2)
        {
            this.DisplacementGradientMode1 = displacementGradientMode1;
            this.DisplacementGradientMode2 = displacementGradientMode2;
            this.StrainTensorMode1 = strainTensorMode1;
            this.StrainTensorMode2 = strainTensorMode2;
            this.StressTensorMode1 = stressTensorMode1;
            this.StressTensorMode2 = stressTensorMode2;
        }
    }
}
