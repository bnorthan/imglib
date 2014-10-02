
package net.imglib2.algorithm.fft2;

import static org.junit.Assert.assertEquals;

import org.jtransforms.fft.FloatFFT_1D;
import org.junit.Test;

import net.imglib2.Cursor;
import net.imglib2.Dimensions;
import net.imglib2.FinalDimensions;
import net.imglib2.FinalInterval;
import net.imglib2.Point;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.FloatArray;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;

/**
 * Test 3rd party fft implementation(s)
 * 
 * @author bnorthan
 */
public class testFFT {

	/**
	 * Tests that the JTransforms api is working within the FFTMethods API
	 */
	@Test
	public void testJTransformFFT1D() {

		int size = 129;

		size = NextSmoothNumber.nextSmooth(size);

		long[] dims = new long[] { size };
		Dimensions dimensions = new FinalDimensions(new long[] { size });

		final float[] array = new float[size];
		final float[] copy = new float[size];

		// place a couple of implules in the array
		array[size / 4] = 50;
		array[size / 2] = 100;

		// make a copy
		for (int i = 0; i < size; i++) {
			copy[i] = array[i];
		}

		// create and perform forward and inverse fft on the data
		final FloatFFT_1D fft = new FloatFFT_1D(size);

		fft.realForward(array);
		fft.realInverse(array, true);

		// assert that we get the original signal back (within an error delta)
		for (int i = 0; i < size; i++) {
			assertEquals(copy[i], array[i], 0.0001);
		}

		// now use the FFTMethods api

		long[] paddedDimensions = new long[1];
		long[] fftDimensions = new long[1];

		// first get the transform dimensions
		FFTMethodsJTransform.dimensionsRealToComplexSmall(dimensions,
			paddedDimensions, fftDimensions);

		// create images for input, transform and inverse
		Img<FloatType> in = ArrayImgs.floats(copy, dims);;

		Img<ComplexFloatType> transform =
			new ArrayImgFactory<ComplexFloatType>().create(fftDimensions,
				new ComplexFloatType());

		Img<FloatType> inverse =
			new ArrayImgFactory<FloatType>().create(dimensions, new FloatType());

		// perform forward and inverse fft
		FFTMethodsJTransform.realToComplex(in, transform, 0);
		FFTMethodsJTransform.complexToReal(transform, inverse, 0);

		int i = 0;

		Cursor<FloatType> cin = in.cursor();
		Cursor<FloatType> cinverse = inverse.cursor();

		while (cin.hasNext()) {

			cin.fwd();
			cinverse.fwd();

			// assert that the inverse = the input within the error delta
			assertEquals(cin.get().getRealFloat(), cinverse.get().getRealFloat(),
				0.0001);

			// assert that the inverse obtained using FFTMethods api is exactly equal
			// to the inverse obtained from using JTransform directly
			assertEquals(array[i], cinverse.get().getRealFloat(), 0);
			i++;
		}
	}

	/**
	 * A test comparing FFTMethods using mines fft with FFTMethods using
	 * JTransform api
	 */
	@Test
	public void testMinesVsJTransform() {

		Dimensions dimensions = new FinalDimensions(new long[] { 100, 100 });

		long[] paddedDimensionsMines = new long[2];
		long[] fftDimensionsMines = new long[2];

		FFTMethods.dimensionsRealToComplexFast(dimensions, paddedDimensionsMines,
			fftDimensionsMines);

		// create an input with a small sphere at the center
		Img<FloatType> in =
			new ArrayImgFactory<FloatType>().create(paddedDimensionsMines,
				new FloatType());
		placeSphereInCenter(in);

		Img<FloatType> inverseMines =
			new ArrayImgFactory<FloatType>().create(paddedDimensionsMines,
				new FloatType());

		Img<FloatType> inverseJTransform =
			new ArrayImgFactory<FloatType>().create(paddedDimensionsMines,
				new FloatType());

		Img<ComplexFloatType> fftMines =
			new ArrayImgFactory<ComplexFloatType>().create(fftDimensionsMines,
				new ComplexFloatType());

		Img<ComplexFloatType> fftJTransform =
			new ArrayImgFactory<ComplexFloatType>().create(fftDimensionsMines,
				new ComplexFloatType());

		FFTMethods.realToComplex(in, fftMines, 0);
		FFTMethodsJTransform.realToComplex(in, fftJTransform, 0);

		FFTMethods.realToComplex(in, fftMines, 0);

		FFTMethodsJTransform.complexToReal(fftJTransform, inverseMines, 0);
		FFTMethodsJTransform.complexToReal(fftJTransform, inverseJTransform, 0);

		Cursor<FloatType> cin = in.cursor();
		Cursor<FloatType> cmines = inverseMines.cursor();
		Cursor<FloatType> cjtransform = inverseJTransform.cursor();

		float inSum = 0;
		float minesSum = 0;
		float jtransformSum = 0;
		int i = 0;

		while (cin.hasNext()) {
			i++;
			cin.fwd();
			cmines.fwd();
			cjtransform.fwd();

			inSum += cin.get().getRealFloat();
			minesSum += cmines.get().getRealFloat();
			jtransformSum += cjtransform.get().getRealFloat();

			assertEquals(cin.get().getRealFloat(), cmines.get().getRealFloat(),
				0.000001);
			assertEquals(cin.get().getRealFloat(), cjtransform.get().getRealFloat(),
				0.000001);
			assertEquals(cjtransform.get().getRealFloat(), cmines.get()
				.getRealFloat(), 0.0);
		}

		System.out.println("FFT even insum: " + inSum + " minessum: " + minesSum);
	}

	// utility to place a small sphere at the center of the image
	private void placeSphereInCenter(Img<FloatType> img) {

		final Point center = new Point(img.numDimensions());

		for (int d = 0; d < img.numDimensions(); d++)
			center.setPosition(img.dimension(d) / 2, d);

		HyperSphere<FloatType> hyperSphere =
			new HyperSphere<FloatType>(img, center, 2);

		for (final FloatType value : hyperSphere) {
			value.setReal(1);
		}
	}

	int seed;

	private int pseudoRandom() {
		return seed = 3170425 * seed + 132102;
	}

	public ArrayImg<FloatType, FloatArray> generateFloatArrayTestImg(
		final boolean fill, final long... dims)
	{
		final float[] array =
			new float[(int) Intervals.numElements(new FinalInterval(dims))];

		if (fill) {
			seed = 17;
			for (int i = 0; i < array.length; i++) {
				array[i] = (float) pseudoRandom() / (float) Integer.MAX_VALUE;
			}
		}

		return ArrayImgs.floats(array, dims);
	}

}
