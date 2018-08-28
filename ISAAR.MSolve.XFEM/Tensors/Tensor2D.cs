using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Tensors
{
    // TODO: Add methods that mustate the same instance instead of returning new ones.
    public class Tensor2D
    {
        public double XX { get; }
        public double YY { get; }
        public double XY { get; }

        public Tensor2D(double componentXX, double componentYY, double componentXY)
        {
            this.XX = componentXX;
            this.YY = componentYY;
            this.XY = componentXY;
        }

        public Tensor2D Multiply(double scalar)
        {
            return new Tensor2D(scalar * XX, scalar * YY, scalar * XY);
        }

        /// <summary>
        /// Colon multiplication of tensors A, B is defined as: A:B = Aij * Bij
        /// </summary>
        /// <param name="tensor"></param>
        /// <returns></returns>
        public double MultiplyColon(Tensor2D tensor)
        {
            return this.XX * tensor.XX + 2 * (this.XY * tensor.XY) + this.YY * tensor.YY;
        }

        /// <summary>
        /// Rotates the tensor counter clockwise.
        /// </summary>
        /// <param name="angle">The counter-clockwise angle that the current tensor is rotated. Should belong to 
        ///     (-π, π] for uniformity, although that isn't checked.</param>
        /// <returns></returns>
        public Tensor2D Rotate(double angle)
        {
            double cos = Math.Cos(2 * angle);
            double sin = Math.Sin(2 * angle);
            double centre = (XX + YY) / 2.0;
            double radius = (XX - YY) / 2.0;

            return new Tensor2D(centre + radius * cos + XY * sin,
                                centre - radius * cos - XY * sin,
                                -radius * sin + XY * cos);
        }

        public Tuple<double, Tensor2D> ToPrincipalAxes()
        {
            throw new NotImplementedException("Decide on the definition interval of the returned angle");
        }

        public Tuple<double, Tensor2D> ToMaxShearAxes()
        {

            throw new NotImplementedException("Decide on the definition interval of the returned angle");
        }
    }
}
