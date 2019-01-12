namespace ISAAR.MSolve.Solvers.LinearSystems
{
    /// <summary>
    /// Objects implementing this interface will be notifying before <see cref="ILinearSystem_v2.Matrix"/> is modified.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ISystemMatrixObserver
    {
        /// <summary>
        /// It will be called before setting <see cref="ILinearSystem_v2.Matrix"/>.
        /// </summary>
        void HandleMatrixWillBeSet();
    }
}
