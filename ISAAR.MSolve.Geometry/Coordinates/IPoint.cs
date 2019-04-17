namespace ISAAR.MSolve.Geometry.Coordinates
{
	/// <summary>
	/// Point in a 3-dimensional space. It can also represent points in 1-dimensional or 2-dimension spaces, but not all 
    /// coordinates will be used. Immutable.
	/// Authors: Serafeim Bakalakos, Dimitris Tsapetis
	/// </summary>
	public interface IPoint
	{
		/// <summary>
		/// Vector with the coordinates of the point. Length = 3.
		/// </summary>
		double[] Coordinates { get; }

        /// <summary>
        /// The first coordinate of the point.
        /// </summary>
        double X1 { get; }

        /// <summary>
        /// The second coordinate of the point.
        /// </summary>
        double X2 { get; }

        /// <summary>
        /// The third coordinate of the point.
        /// </summary>
        double X3 { get; }
    }
}