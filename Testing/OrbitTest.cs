using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Testing
{
    [TestClass]
    public class OrbitTest
    {
        [TestMethod]
        public void ConstructorTest()
        {
            float a = 1, e = 0, M = 1;
            Orbit orbit = new(a ,e, M);
            Assert.AreEqual(orbit.SemiMajorAxis(), a);
            Assert.AreEqual(orbit.Eccentricity(), e);
            Assert.AreEqual(orbit.StandardGrav(), M);
        }
    }

}
