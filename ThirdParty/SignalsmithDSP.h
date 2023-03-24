#pragma once 
#include <complex>
#include <cmath>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

namespace signalsmith
{
    // Perf.h
    namespace perf {
	/**	@defgroup Performance Performance helpers
		@brief Nothing serious, just some `#defines` and helpers
		
		@{
		@file
	*/
		
	/// *Really* insist that a function/method is inlined (mostly for performance in DEBUG builds)
	#ifndef SIGNALSMITH_INLINE
	#ifdef __GNUC__
	#define SIGNALSMITH_INLINE __attribute__((always_inline)) inline
	#elif defined(__MSVC__)
	#define SIGNALSMITH_INLINE __forceinline inline
	#else
	#define SIGNALSMITH_INLINE inline
	#endif
	#endif

	/** @brief Complex-multiplication (with optional conjugate second-arg), without handling NaN/Infinity
		The `std::complex` multiplication has edge-cases around NaNs which slow things down and prevent auto-vectorisation.
	*/
	template <bool conjugateSecond=false, typename V>
	SIGNALSMITH_INLINE static std::complex<V> mul(const std::complex<V> &a, const std::complex<V> &b) {
		return conjugateSecond ? std::complex<V>{
			b.real()*a.real() + b.imag()*a.imag(),
			b.real()*a.imag() - b.imag()*a.real()
		} : std::complex<V>{
			a.real()*b.real() - a.imag()*b.imag(),
			a.real()*b.imag() + a.imag()*b.real()
		};
	}
}
    // fft.h
    namespace fft {
	/**	@defgroup FFT FFT (complex and real)
		@brief Fourier transforms (complex and real)

		@{
		@file
	*/

	namespace _fft_impl {

		template <typename V>
		SIGNALSMITH_INLINE V complexReal(const std::complex<V> &c) {
			return ((V*)(&c))[0];
		}
		template <typename V>
		SIGNALSMITH_INLINE V complexImag(const std::complex<V> &c) {
			return ((V*)(&c))[1];
		}

		// Complex multiplication has edge-cases around Inf/NaN - handling those properly makes std::complex non-inlineable, so we use our own
		template <bool conjugateSecond, typename V>
		SIGNALSMITH_INLINE std::complex<V> complexMul(const std::complex<V> &a, const std::complex<V> &b) {
			V aReal = complexReal(a), aImag = complexImag(a);
			V bReal = complexReal(b), bImag = complexImag(b);
			return conjugateSecond ? std::complex<V>{
				bReal*aReal + bImag*aImag,
				bReal*aImag - bImag*aReal
			} : std::complex<V>{
				aReal*bReal - aImag*bImag,
				aReal*bImag + aImag*bReal
			};
		}

		template<bool flipped, typename V>
		SIGNALSMITH_INLINE std::complex<V> complexAddI(const std::complex<V> &a, const std::complex<V> &b) {
			V aReal = complexReal(a), aImag = complexImag(a);
			V bReal = complexReal(b), bImag = complexImag(b);
			return flipped ? std::complex<V>{
				aReal + bImag,
				aImag - bReal
			} : std::complex<V>{
				aReal - bImag,
				aImag + bReal
			};
		}

		// Use SFINAE to get an iterator from std::begin(), if supported - otherwise assume the value itself is an iterator
		template<typename T, typename=void>
		struct GetIterator {
			static T get(const T &t) {
				return t;
			}
		};
		template<typename T>
		struct GetIterator<T, decltype((void)std::begin(std::declval<T>()))> {
			static auto get(const T &t) -> decltype(std::begin(t)) {
				return std::begin(t);
			}
		};
	}

	/** Floating-point FFT implementation.
	It is fast for 2^a * 3^b.
	Here are the peak and RMS errors for `float`/`double` computation:
	\diagram{fft-errors.svg Simulated errors for pure-tone harmonic inputs\, compared to a theoretical upper bound from "Roundoff error analysis of the fast Fourier transform" (G. Ramos, 1971)}
	*/
	template<typename V=double>
	class FFT {
		using complex = std::complex<V>;
		size_t _size;
		std::vector<complex> workingVector;
		
		enum class StepType {
			generic, step2, step3, step4
		};
		struct Step {
			StepType type;
			size_t factor;
			size_t startIndex;
			size_t innerRepeats;
			size_t outerRepeats;
			size_t twiddleIndex;
		};
		std::vector<size_t> factors;
		std::vector<Step> plan;
		std::vector<complex> twiddleVector;
		
		struct PermutationPair {size_t from, to;};
		std::vector<PermutationPair> permutation;
		
		void addPlanSteps(size_t factorIndex, size_t start, size_t length, size_t repeats) {
			if (factorIndex >= factors.size()) return;
			
			size_t factor = factors[factorIndex];
			if (factorIndex + 1 < factors.size()) {
				if (factors[factorIndex] == 2 && factors[factorIndex + 1] == 2) {
					++factorIndex;
					factor = 4;
				}
			}

			size_t subLength = length/factor;
			Step mainStep{StepType::generic, factor, start, subLength, repeats, twiddleVector.size()};

			if (factor == 2) mainStep.type = StepType::step2;
			if (factor == 3) mainStep.type = StepType::step3;
			if (factor == 4) mainStep.type = StepType::step4;

			// Twiddles
			bool foundStep = false;
			for (const Step &existingStep : plan) {
				if (existingStep.factor == mainStep.factor && existingStep.innerRepeats == mainStep.innerRepeats) {
					foundStep = true;
					mainStep.twiddleIndex = existingStep.twiddleIndex;
					break;
				}
			}
			if (!foundStep) {
				for (size_t i = 0; i < subLength; ++i) {
					for (size_t f = 0; f < factor; ++f) {
						double phase = 2*M_PI*i*f/length;
						complex twiddle = {V(std::cos(phase)), V(-std::sin(phase))};
						twiddleVector.push_back(twiddle);
					}
				}
			}

			if (repeats == 1 && sizeof(complex)*subLength > 65536) {
				for (size_t i = 0; i < factor; ++i) {
					addPlanSteps(factorIndex + 1, start + i*subLength, subLength, 1);
				}
			} else {
				addPlanSteps(factorIndex + 1, start, subLength, repeats*factor);
			}
			plan.push_back(mainStep);
		}
		void setPlan() {
			factors.resize(0);
			size_t size = _size, factor = 2;
			while (size > 1) {
				if (size%factor == 0) {
					factors.push_back(factor);
					size /= factor;
				} else if (factor > sqrt(size)) {
					factor = size;
				} else {
					++factor;
				}
			}

			plan.resize(0);
			twiddleVector.resize(0);
			addPlanSteps(0, 0, _size, 1);
			
			permutation.resize(0);
			permutation.push_back(PermutationPair{0, 0});
			size_t indexLow = 0, indexHigh = factors.size();
			size_t inputStepLow = _size, outputStepLow = 1;
			size_t inputStepHigh = 1, outputStepHigh = _size;
			while (outputStepLow*inputStepHigh < _size) {
				size_t f, inputStep, outputStep;
				if (outputStepLow <= inputStepHigh) {
					f = factors[indexLow++];
					inputStep = (inputStepLow /= f);
					outputStep = outputStepLow;
					outputStepLow *= f;
				} else {
					f = factors[--indexHigh];
					inputStep = inputStepHigh;
					inputStepHigh *= f;
					outputStep = (outputStepHigh /= f);
				}
				size_t oldSize = permutation.size();
				for (size_t i = 1; i < f; ++i) {
					for (size_t j = 0; j < oldSize; ++j) {
						PermutationPair pair = permutation[j];
						pair.from += i*inputStep;
						pair.to += i*outputStep;
						permutation.push_back(pair);
					}
				}
			}
		}

		template<bool inverse, typename RandomAccessIterator>
		void fftStepGeneric(RandomAccessIterator &&origData, const Step &step) {
			complex *working = workingVector.data();
			const size_t stride = step.innerRepeats;

			for (size_t outerRepeat = 0; outerRepeat < step.outerRepeats; ++outerRepeat) {
				RandomAccessIterator data = origData;
				
				const complex *twiddles = twiddleVector.data() + step.twiddleIndex;
				const size_t factor = step.factor;
				for (size_t repeat = 0; repeat < step.innerRepeats; ++repeat) {
					for (size_t i = 0; i < step.factor; ++i) {
						working[i] = _fft_impl::complexMul<inverse>(data[i*stride], twiddles[i]);
					}
					for (size_t f = 0; f < factor; ++f) {
						complex sum = working[0];
						for (size_t i = 1; i < factor; ++i) {
							double phase = 2*M_PI*f*i/factor;
							complex twiddle = {V(std::cos(phase)), V(-std::sin(phase))};
							sum += _fft_impl::complexMul<inverse>(working[i], twiddle);
						}
						data[f*stride] = sum;
					}
					++data;
					twiddles += factor;
				}
				origData += step.factor*step.innerRepeats;
			}
		}

		template<bool inverse, typename RandomAccessIterator>
		SIGNALSMITH_INLINE void fftStep2(RandomAccessIterator &&origData, const Step &step) {
			const size_t stride = step.innerRepeats;
			const complex *origTwiddles = twiddleVector.data() + step.twiddleIndex;
			for (size_t outerRepeat = 0; outerRepeat < step.outerRepeats; ++outerRepeat) {
				const complex* twiddles = origTwiddles;
				for (RandomAccessIterator data = origData; data < origData + stride; ++data) {
					complex A = data[0];
					complex B = _fft_impl::complexMul<inverse>(data[stride], twiddles[1]);
					
					data[0] = A + B;
					data[stride] = A - B;
					twiddles += 2;
				}
				origData += 2*stride;
			}
		}

		template<bool inverse, typename RandomAccessIterator>
		SIGNALSMITH_INLINE void fftStep3(RandomAccessIterator &&origData, const Step &step) {
			constexpr complex factor3 = {-0.5, inverse ? 0.8660254037844386 : -0.8660254037844386};
			const size_t stride = step.innerRepeats;
			const complex *origTwiddles = twiddleVector.data() + step.twiddleIndex;
			
			for (size_t outerRepeat = 0; outerRepeat < step.outerRepeats; ++outerRepeat) {
				const complex* twiddles = origTwiddles;
				for (RandomAccessIterator data = origData; data < origData + stride; ++data) {
					complex A = data[0];
					complex B = _fft_impl::complexMul<inverse>(data[stride], twiddles[1]);
					complex C = _fft_impl::complexMul<inverse>(data[stride*2], twiddles[2]);
					
					complex realSum = A + (B + C)*factor3.real();
					complex imagSum = (B - C)*factor3.imag();

					data[0] = A + B + C;
					data[stride] = _fft_impl::complexAddI<false>(realSum, imagSum);
					data[stride*2] = _fft_impl::complexAddI<true>(realSum, imagSum);

					twiddles += 3;
				}
				origData += 3*stride;
			}
		}

		template<bool inverse, typename RandomAccessIterator>
		SIGNALSMITH_INLINE void fftStep4(RandomAccessIterator &&origData, const Step &step) {
			const size_t stride = step.innerRepeats;
			const complex *origTwiddles = twiddleVector.data() + step.twiddleIndex;
			
			for (size_t outerRepeat = 0; outerRepeat < step.outerRepeats; ++outerRepeat) {
				const complex* twiddles = origTwiddles;
				for (RandomAccessIterator data = origData; data < origData + stride; ++data) {
					complex A = data[0];
					complex C = _fft_impl::complexMul<inverse>(data[stride], twiddles[2]);
					complex B = _fft_impl::complexMul<inverse>(data[stride*2], twiddles[1]);
					complex D = _fft_impl::complexMul<inverse>(data[stride*3], twiddles[3]);

					complex sumAC = A + C, sumBD = B + D;
					complex diffAC = A - C, diffBD = B - D;

					data[0] = sumAC + sumBD;
					data[stride] = _fft_impl::complexAddI<!inverse>(diffAC, diffBD);
					data[stride*2] = sumAC - sumBD;
					data[stride*3] = _fft_impl::complexAddI<inverse>(diffAC, diffBD);

					twiddles += 4;
				}
				origData += 4*stride;
			}
		}
		
		template<typename InputIterator, typename OutputIterator>
		void permute(InputIterator input, OutputIterator data) {
			for (auto pair : permutation) {
				data[pair.from] = input[pair.to];
			}
		}

		template<bool inverse, typename InputIterator, typename OutputIterator>
		void run(InputIterator &&input, OutputIterator &&data) {
			permute(input, data);
			
			for (const Step &step : plan) {
				switch (step.type) {
					case StepType::generic:
						fftStepGeneric<inverse>(data + step.startIndex, step);
						break;
					case StepType::step2:
						fftStep2<inverse>(data + step.startIndex, step);
						break;
					case StepType::step3:
						fftStep3<inverse>(data + step.startIndex, step);
						break;
					case StepType::step4:
						fftStep4<inverse>(data + step.startIndex, step);
						break;
				}
			}
		}

		static bool validSize(size_t size) {
			constexpr static bool filter[32] = {
				1, 1, 1, 1, 1, 0, 1, 0, 1, 1, // 0-9
				0, 0, 1, 0, 0, 0, 1, 0, 1, 0, // 10-19
				0, 0, 0, 0, 1, 0, 0, 0, 0, 0, // 20-29
				0, 0
			};
			return filter[size];
		}
	public:
		static size_t fastSizeAbove(size_t size) {
			size_t power2 = 1;
			while (size >= 32) {
				size = (size - 1)/2 + 1;
				power2 *= 2;
			}
			while (size < 32 && !validSize(size)) {
				++size;
			}
			return power2*size;
		}
		static size_t fastSizeBelow(size_t size) {
			size_t power2 = 1;
			while (size >= 32) {
				size /= 2;
				power2 *= 2;
			}
			while (size > 1 && !validSize(size)) {
				--size;
			}
			return power2*size;
		}

		FFT(size_t size, int fastDirection=0) : _size(0) {
			if (fastDirection > 0) size = fastSizeAbove(size);
			if (fastDirection < 0) size = fastSizeBelow(size);
			this->setSize(size);
		}

		size_t setSize(size_t size) {
			if (size != _size) {
				_size = size;
				workingVector.resize(size);
				setPlan();
			}
			return _size;
		}
		size_t setFastSizeAbove(size_t size) {
			return setSize(fastSizeAbove(size));
		}
		size_t setFastSizeBelow(size_t size) {
			return setSize(fastSizeBelow(size));
		}
		const size_t & size() const {
			return _size;
		}

		template<typename InputIterator, typename OutputIterator>
		void fft(InputIterator &&input, OutputIterator &&output) {
			auto inputIter = _fft_impl::GetIterator<InputIterator>::get(input);
			auto outputIter = _fft_impl::GetIterator<OutputIterator>::get(output);
			return run<false>(inputIter, outputIter);
		}

		template<typename InputIterator, typename OutputIterator>
		void ifft(InputIterator &&input, OutputIterator &&output) {
			auto inputIter = _fft_impl::GetIterator<InputIterator>::get(input);
			auto outputIter = _fft_impl::GetIterator<OutputIterator>::get(output);
			return run<true>(inputIter, outputIter);
		}
	};

	struct FFTOptions {
		static constexpr int halfFreqShift = 1;
	};

	template<typename V, int optionFlags=0>
	class RealFFT {
		static constexpr bool modified = (optionFlags&FFTOptions::halfFreqShift);

		using complex = std::complex<V>;
		std::vector<complex> complexBuffer1, complexBuffer2;
		std::vector<complex> twiddlesMinusI;
		std::vector<complex> modifiedRotations;
		FFT<V> complexFft;
	public:
		static size_t fastSizeAbove(size_t size) {
			return FFT<V>::fastSizeAbove((size + 1)/2)*2;
		}
		static size_t fastSizeBelow(size_t size) {
			return FFT<V>::fastSizeBelow(size/2)*2;
		}

		RealFFT(size_t size=0, int fastDirection=0) : complexFft(0) {
			if (fastDirection > 0) size = fastSizeAbove(size);
			if (fastDirection < 0) size = fastSizeBelow(size);
			this->setSize(std::max<size_t>(size, 2));
		}

		size_t setSize(size_t size) {
			complexBuffer1.resize(size/2);
			complexBuffer2.resize(size/2);

			size_t hhSize = size/4 + 1;
			twiddlesMinusI.resize(hhSize);
			for (size_t i = 0; i < hhSize; ++i) {
				V rotPhase = -2*M_PI*(modified ? i + 0.5 : i)/size;
				twiddlesMinusI[i] = {std::sin(rotPhase), -std::cos(rotPhase)};
			}
			if (modified) {
				modifiedRotations.resize(size/2);
				for (size_t i = 0; i < size/2; ++i) {
					V rotPhase = -2*M_PI*i/size;
					modifiedRotations[i] = {std::cos(rotPhase), std::sin(rotPhase)};
				}
			}
			
			return complexFft.setSize(size/2);
		}
		size_t setFastSizeAbove(size_t size) {
			return setSize(fastSizeAbove(size));
		}
		size_t setFastSizeBelow(size_t size) {
			return setSize(fastSizeBelow(size));
		}
		size_t size() const {
			return complexFft.size()*2;
		}

		template<typename InputIterator, typename OutputIterator>
		void fft(InputIterator &&input, OutputIterator &&output) {
			size_t hSize = complexFft.size();
			for (size_t i = 0; i < hSize; ++i) {
				if (modified) {
					complexBuffer1[i] = _fft_impl::complexMul<false>({input[2*i], input[2*i + 1]}, modifiedRotations[i]);
				} else {
					complexBuffer1[i] = {input[2*i], input[2*i + 1]};
				}
			}
			
			complexFft.fft(complexBuffer1.data(), complexBuffer2.data());
			
			if (!modified) output[0] = {
				complexBuffer2[0].real() + complexBuffer2[0].imag(),
				complexBuffer2[0].real() - complexBuffer2[0].imag()
			};
			for (size_t i = modified ? 0 : 1; i <= hSize/2; ++i) {
				size_t conjI = modified ? (hSize  - 1 - i) : (hSize - i);
				
				complex odd = (complexBuffer2[i] + conj(complexBuffer2[conjI]))*(V)0.5;
				complex evenI = (complexBuffer2[i] - conj(complexBuffer2[conjI]))*(V)0.5;
				complex evenRotMinusI = _fft_impl::complexMul<false>(evenI, twiddlesMinusI[i]);

				output[i] = odd + evenRotMinusI;
				output[conjI] = conj(odd - evenRotMinusI);
			}
		}

		template<typename InputIterator, typename OutputIterator>
		void ifft(InputIterator &&input, OutputIterator &&output) {
			size_t hSize = complexFft.size();
			if (!modified) complexBuffer1[0] = {
				input[0].real() + input[0].imag(),
				input[0].real() - input[0].imag()
			};
			for (size_t i = modified ? 0 : 1; i <= hSize/2; ++i) {
				size_t conjI = modified ? (hSize  - 1 - i) : (hSize - i);
				complex v = input[i], v2 = input[conjI];

				complex odd = v + conj(v2);
				complex evenRotMinusI = v - conj(v2);
				complex evenI = _fft_impl::complexMul<true>(evenRotMinusI, twiddlesMinusI[i]);
				
				complexBuffer1[i] = odd + evenI;
				complexBuffer1[conjI] = conj(odd - evenI);
			}
			
			complexFft.ifft(complexBuffer1.data(), complexBuffer2.data());
			
			for (size_t i = 0; i < hSize; ++i) {
				complex v = complexBuffer2[i];
				if (modified) v = _fft_impl::complexMul<true>(v, modifiedRotations[i]);
				output[2*i] = v.real();
				output[2*i + 1] = v.imag();
			}
		}
	};

	template<typename V>
	struct ModifiedRealFFT : public RealFFT<V, FFTOptions::halfFreqShift> {
		using RealFFT<V, FFTOptions::halfFreqShift>::RealFFT;
	};

/// @}
}
    // windows.h
    namespace windows {
	/**	@defgroup Windows Window functions
		@brief Windows for spectral analysis
		
		These are generally double-precision, because they are mostly calculated during setup/reconfiguring, not real-time code.
		
		@{
		@file
	*/
	
	/** @brief The Kaiser window (almost) maximises the energy in the main-lobe compared to the side-lobes.
		
		Kaiser windows can be constructing using the shape-parameter (beta) or using the static `with???()` methods.*/
	class Kaiser {
		// I_0(x)=\sum_{k=0}^{N}\frac{x^{2k}}{(k!)^2\cdot4^k}
		inline static double bessel0(double x) {
			const double significanceLimit = 1e-4;
			double result = 0;
			double term = 1;
			double m = 0;
			while (term > significanceLimit) {
				result += term;
				++m;
				term *= (x*x)/(4*m*m);
			}

			return result;
		}
		double beta;
		double invB0;
		
		static double heuristicBandwidth(double bandwidth) {
			// Good peaks
			//return bandwidth + 8/((bandwidth + 3)*(bandwidth + 3));
			// Good average
			//return bandwidth + 14/((bandwidth + 2.5)*(bandwidth + 2.5));
			// Compromise
			return bandwidth + 8/((bandwidth + 3)*(bandwidth + 3)) + 0.25*std::max(3 - bandwidth, 0.0);
		}
	public:
		/// Set up a Kaiser window with a given shape.  `beta` is `pi*alpha` (since there is ambiguity about shape parameters)
		Kaiser(double beta) : beta(beta), invB0(1/bessel0(beta)) {}

		/// @name Bandwidth methods
		/// @{
		static Kaiser withBandwidth(double bandwidth, bool heuristicOptimal=false) {
			return Kaiser(bandwidthToBeta(bandwidth, heuristicOptimal));
		}

		/** Returns the Kaiser shape where the main lobe has the specified bandwidth (as a factor of 1/window-length).
		\diagram{kaiser-windows.svg,You can see that the main lobe matches the specified bandwidth.}
		If `heuristicOptimal` is enabled, the main lobe width is _slightly_ wider, improving both the peak and total energy - see `bandwidthToEnergyDb()` and `bandwidthToPeakDb()`.
		\diagram{kaiser-windows-heuristic.svg, The main lobe extends to ±bandwidth/2.} */
		static double bandwidthToBeta(double bandwidth, bool heuristicOptimal=false) {
			if (heuristicOptimal) { // Heuristic based on numerical search
				bandwidth = heuristicBandwidth(bandwidth);
			}
			bandwidth = std::max(bandwidth, 2.0);
			double alpha = std::sqrt(bandwidth*bandwidth*0.25 - 1);
			return alpha*M_PI;
		}
		
		static double betaToBandwidth(double beta) {
			double alpha = beta*(1.0/M_PI);
			return 2*std::sqrt(alpha*alpha + 1);
		}
		/// @}

		/// @name Performance methods
		/// @{
		/** @brief Total energy ratio (in dB) between side-lobes and the main lobe.
			\diagram{windows-kaiser-sidelobe-energy.svg,Measured main/side lobe energy ratio.  You can see that the heuristic improves performance for all bandwidth values.}
			This function uses an approximation which is accurate to ±0.5dB for 2 ⩽ bandwidth ≤ 10, or 1 ⩽ bandwidth ≤ 10 when `heuristicOptimal`is enabled.
		*/
		static double bandwidthToEnergyDb(double bandwidth, bool heuristicOptimal=false) {
			// Horrible heuristic fits
			if (heuristicOptimal) {
				if (bandwidth < 3) bandwidth += (3 - bandwidth)*0.5;
				return 12.9 + -3/(bandwidth + 0.4) - 13.4*bandwidth + (bandwidth < 3)*-9.6*(bandwidth - 3);
			}
			return 10.5 + 15/(bandwidth + 0.4) - 13.25*bandwidth + (bandwidth < 2)*13*(bandwidth - 2);
		}
		static double energyDbToBandwidth(double energyDb, bool heuristicOptimal=false) {
			double bw = 1;
			while (bw < 20 && bandwidthToEnergyDb(bw, heuristicOptimal) > energyDb) {
				bw *= 2;
			}
			double step = bw/2;
			while (step > 0.0001) {
				if (bandwidthToEnergyDb(bw, heuristicOptimal) > energyDb) {
					bw += step;
				} else {
					bw -= step;
				}
				step *= 0.5;
			}
			return bw;
		}
		/** @brief Peak ratio (in dB) between side-lobes and the main lobe.
			\diagram{windows-acg-sidelobe-peaks.svg,Measured main/side lobe peak ratio.  You can see that the heuristic improves performance, except in the bandwidth range 1-2 where peak ratio was sacrificed to improve total energy ratio.}
			This function uses an approximation which is accurate to ±0.5dB for 2 ⩽ bandwidth ≤ 9, or 0.5 ⩽ bandwidth ≤ 9 when `heuristicOptimal`is enabled.
		*/
		static double bandwidthToPeakDb(double bandwidth, bool heuristicOptimal=false) {
			// Horrible heuristic fits
			if (heuristicOptimal) {
				return 14.2 - 20/(bandwidth + 1) - 13*bandwidth + (bandwidth < 3)*-6*(bandwidth - 3) + (bandwidth < 2.25)*5.8*(bandwidth - 2.25);
			}
			return 10 + 8/(bandwidth + 2) - 12.75*bandwidth + (bandwidth < 2)*4*(bandwidth - 2);
		}
		static double peakDbToBandwidth(double peakDb, bool heuristicOptimal=false) {
			double bw = 1;
			while (bw < 20 && bandwidthToPeakDb(bw, heuristicOptimal) > peakDb) {
				bw *= 2;
			}
			double step = bw/2;
			while (step > 0.0001) {
				if (bandwidthToPeakDb(bw, heuristicOptimal) > peakDb) {
					bw += step;
				} else {
					bw -= step;
				}
				step *= 0.5;
			}
			return bw;
		}
		/** @} */

		/** Equivalent noise bandwidth (ENBW), a measure of frequency resolution.
			\diagram{windows-kaiser-enbw.svg,Measured ENBW\, with and without the heuristic bandwidth adjustment.}
			This approximation is accurate to ±0.05 up to a bandwidth of 22.
		*/
		static double bandwidthToEnbw(double bandwidth, bool heuristicOptimal=false) {
			if (heuristicOptimal) bandwidth = heuristicBandwidth(bandwidth);
			double b2 = std::max<double>(bandwidth - 2, 0);
			return 1 + b2*(0.2 + b2*(-0.005 + b2*(-0.000005 + b2*0.0000022)));
		}

		/// Return the window's value for position in the range [0, 1]
		double operator ()(double unit) {
			double r = 2*unit - 1;
			double arg = std::sqrt(1 - r*r);
			return bessel0(beta*arg)*invB0;
		}
	
		/// Fills an arbitrary container with a Kaiser window
		template<typename Data>
		void fill(Data &data, int size) const {
			double invSize = 1.0/size;
			for (int i = 0; i < size; ++i) {
				double r = (2*i + 1)*invSize - 1;
				double arg = std::sqrt(1 - r*r);
				data[i] = bessel0(beta*arg)*invB0;
			}
		}
	};

	/** @brief The Approximate Confined Gaussian window is (almost) optimal
		
		ACG windows can be constructing using the shape-parameter (sigma) or using the static `with???()` methods.*/
	class ApproximateConfinedGaussian {
		double gaussianFactor;
		
		double gaussian(double x) const {
			return std::exp(-x*x*gaussianFactor);
		}
	public:
		/// Heuristic map from bandwidth to the appropriately-optimal sigma
		static double bandwidthToSigma(double bandwidth) {
			return 0.3/std::sqrt(bandwidth);
		}
		static ApproximateConfinedGaussian withBandwidth(double bandwidth) {
			return ApproximateConfinedGaussian(bandwidthToSigma(bandwidth));
		}

		ApproximateConfinedGaussian(double sigma) : gaussianFactor(0.0625/(sigma*sigma)) {}
	
		/// Fills an arbitrary container
		template<typename Data>
		void fill(Data &data, int size) const {
			double invSize = 1.0/size;
			double offsetScale = gaussian(1)/(gaussian(3) + gaussian(-1));
			double norm = 1/(gaussian(0) - 2*offsetScale*(gaussian(2)));
			for (int i = 0; i < size; ++i) {
				double r = (2*i + 1)*invSize - 1;
				data[i] = norm*(gaussian(r) - offsetScale*(gaussian(r - 2) + gaussian(r + 2)));
			}
		}
	};

	/** Forces STFT perfect-reconstruction (WOLA) on an existing window, for a given STFT interval.
	For example, here are perfect-reconstruction versions of the approximately-optimal @ref Kaiser windows:
	\diagram{kaiser-windows-heuristic-pr.svg,Note the lower overall energy\, and the pointy top for 2x bandwidth. Spectral performance is about the same\, though.}
	*/
	template<typename Data>
	void forcePerfectReconstruction(Data &data, int windowLength, int interval) {
		for (int i = 0; i < interval; ++i) {
			double sum2 = 0;
			for (int index = i; index < windowLength; index += interval) {
				sum2 += data[index]*data[index];
			}
			double factor = 1/std::sqrt(sum2);
			for (int index = i; index < windowLength; index += interval) {
				data[index] *= factor;
			}
		}
	}

/** @} */
}
    // delay.h
    namespace delay {
	/**	@defgroup PitchDelayProcessor PitchDelayProcessor utilities
		@brief Standalone templated classes for delays
		
		You can set up a `Buffer` or `MultiBuffer`, and get interpolated samples using a `Reader` (separately on each channel in the multi-channel case) - or you can use `PitchDelayProcessor`/`MultiDelay` which include their own buffers.

		Interpolation quality is chosen using a template class, from @ref Interpolators.

		@{
		@file
	*/

	/** @brief Single-channel delay buffer
 
		Access is used with `buffer[]`, relative to the internal read/write position ("head").  This head is moved using `++buffer` (or `buffer += n`), such that `buffer[1] == (buffer + 1)[0]` in a similar way iterators/pointers.
		
		Operations like `buffer - 10` or `buffer++` return a View, which holds a fixed position in the buffer (based on the read/write position at the time).
		
		The capacity includes both positive and negative indices.  For example, a capacity of 100 would support using any of the ranges:
		
		* `buffer[-99]` to buffer[0]`
		* `buffer[-50]` to buffer[49]`
		* `buffer[0]` to buffer[99]`

		Although buffers are usually used with historical samples accessed using negative indices e.g. `buffer[-10]`, you could equally use it flipped around (moving the head backwards through the buffer using `--buffer`).
	*/
	template<typename Sample>
	class Buffer {
		unsigned bufferIndex;
		unsigned bufferMask;
		std::vector<Sample> buffer;
	public:
		Buffer(int minCapacity=0) {
			resize(minCapacity);
		}
		// We shouldn't accidentally copy a delay buffer
		Buffer(const Buffer &other) = delete;
		Buffer & operator =(const Buffer &other) = delete;
		// But moving one is fine
		Buffer(Buffer &&other) = default;
		Buffer & operator =(Buffer &&other) = default;

		void resize(int minCapacity, Sample value=Sample()) {
			int bufferLength = 1;
			while (bufferLength < minCapacity) bufferLength *= 2;
			buffer.assign(bufferLength, value);
			bufferMask = unsigned(bufferLength - 1);
			bufferIndex = 0;
		}
		void reset(Sample value=Sample()) {
			buffer.assign(buffer.size(), value);
		}

		/// Holds a view for a particular position in the buffer
		template<bool isConst>
		class View {
			using CBuffer = typename std::conditional<isConst, const Buffer, Buffer>::type;
			using CSample = typename std::conditional<isConst, const Sample, Sample>::type;
			CBuffer *buffer = nullptr;
			unsigned bufferIndex = 0;
		public:
			View(CBuffer &buffer, int offset=0) : buffer(&buffer), bufferIndex(buffer.bufferIndex + (unsigned)offset) {}
			View(const View &other, int offset=0) : buffer(other.buffer), bufferIndex(other.bufferIndex + (unsigned)offset) {}
			View & operator =(const View &other) {
				buffer = other.buffer;
				bufferIndex = other.bufferIndex;
				return *this;
			}
			
			CSample & operator[](int offset) {
				return buffer->buffer[(bufferIndex + (unsigned)offset)&buffer->bufferMask];
			}
			const Sample & operator[](int offset) const {
				return buffer->buffer[(bufferIndex + (unsigned)offset)&buffer->bufferMask];
			}

			/// Write data into the buffer
			template<typename Data>
			void write(Data &&data, int length) {
				for (int i = 0; i < length; ++i) {
					(*this)[i] = data[i];
				}
			}
			/// Read data out from the buffer
			template<typename Data>
			void read(int length, Data &&data) const {
				for (int i = 0; i < length; ++i) {
					data[i] = (*this)[i];
				}
			}

			View operator +(int offset) const {
				return View(*this, offset);
			}
			View operator -(int offset) const {
				return View(*this, -offset);
			}
		};
		using MutableView = View<false>;
		using ConstView = View<true>;
		
		MutableView view(int offset=0) {
			return MutableView(*this, offset);
		}
		ConstView view(int offset=0) const {
			return ConstView(*this, offset);
		}
		ConstView constView(int offset=0) const {
			return ConstView(*this, offset);
		}

		Sample & operator[](int offset) {
			return buffer[(bufferIndex + (unsigned)offset)&bufferMask];
		}
		const Sample & operator[](int offset) const {
			return buffer[(bufferIndex + (unsigned)offset)&bufferMask];
		}

		/// Write data into the buffer
		template<typename Data>
		void write(Data &&data, int length) {
			for (int i = 0; i < length; ++i) {
				(*this)[i] = data[i];
			}
		}
		/// Read data out from the buffer
		template<typename Data>
		void read(int length, Data &&data) const {
			for (int i = 0; i < length; ++i) {
				data[i] = (*this)[i];
			}
		}
		
		Buffer & operator ++() {
			++bufferIndex;
			return *this;
		}
		Buffer & operator +=(int i) {
			bufferIndex += (unsigned)i;
			return *this;
		}
		Buffer & operator --() {
			--bufferIndex;
			return *this;
		}
		Buffer & operator -=(int i) {
			bufferIndex -= (unsigned)i;
			return *this;
		}

		MutableView operator ++(int) {
			MutableView view(*this);
			++bufferIndex;
			return view;
		}
		MutableView operator +(int i) {
			return MutableView(*this, i);
		}
		ConstView operator +(int i) const {
			return ConstView(*this, i);
		}
		MutableView operator --(int) {
			MutableView view(*this);
			--bufferIndex;
			return view;
		}
		MutableView operator -(int i) {
			return MutableView(*this, -i);
		}
		ConstView operator -(int i) const {
			return ConstView(*this, -i);
		}
	};

	/** @brief Multi-channel delay buffer

		This behaves similarly to the single-channel `Buffer`, with the following differences:
		
		* `buffer[c]` returns a view for a single channel, which behaves like the single-channel `Buffer::View`.
		* The constructor and `.resize()` take an additional first `channel` argument.
	*/
	template<typename Sample>
	class MultiBuffer {
		int channels, stride;
		Buffer<Sample> buffer;
	public:
		using ConstChannel = typename Buffer<Sample>::ConstView;
		using MutableChannel = typename Buffer<Sample>::MutableView;

		MultiBuffer(int channels=0, int capacity=0) : channels(channels), stride(capacity), buffer(channels*capacity) {}

		void resize(int nChannels, int capacity, Sample value=Sample()) {
			channels = nChannels;
			stride = capacity;
			buffer.resize(channels*capacity, value);
		}
		void reset(Sample value=Sample()) {
			buffer.reset(value);
		}

		/// A reference-like multi-channel result for a particular sample index
		template<bool isConst>
		class Stride {
			using CChannel = typename std::conditional<isConst, ConstChannel, MutableChannel>::type;
			using CSample = typename std::conditional<isConst, const Sample, Sample>::type;
			CChannel view;
			int channels, stride;
		public:
			Stride(CChannel view, int channels, int stride) : view(view), channels(channels), stride(stride) {}
			Stride(const Stride &other) : view(other.view), channels(other.channels), stride(other.stride) {}
			
			CSample & operator[](int channel) {
				return view[channel*stride];
			}
			const Sample & operator[](int channel) const {
				return view[channel*stride];
			}

			/// Reads from the buffer into a multi-channel result
			template<class Data>
			void get(Data &&result) const {
				for (int c = 0; c < channels; ++c) {
					result[c] = view[c*stride];
				}
			}
			/// Writes from multi-channel data into the buffer
			template<class Data>
			void set(Data &&data) {
				for (int c = 0; c < channels; ++c) {
					view[c*stride] = data[c];
				}
			}
			template<class Data>
			Stride & operator =(const Data &data) {
				set(data);
				return *this;
			}
			Stride & operator =(const Stride &data) {
				set(data);
				return *this;
			}
		};
		
		Stride<false> at(int offset) {
			return {buffer.view(offset), channels, stride};
		}
		Stride<true> at(int offset) const {
			return {buffer.view(offset), channels, stride};
		}

		/// Holds a particular position in the buffer
		template<bool isConst>
		class View {
			using CChannel = typename std::conditional<isConst, ConstChannel, MutableChannel>::type;
			CChannel view;
			int channels, stride;
		public:
			View(CChannel view, int channels, int stride) : view(view), channels(channels), stride(stride) {}
			
			CChannel operator[](int channel) {
				return view + channel*stride;
			}
			ConstChannel operator[](int channel) const {
				return view + channel*stride;
			}

			Stride<isConst> at(int offset) {
				return {view + offset, channels, stride};
			}
			Stride<true> at(int offset) const {
				return {view + offset, channels, stride};
			}
		};
		using MutableView = View<false>;
		using ConstView = View<true>;

		MutableView view(int offset=0) {
			return MutableView(buffer.view(offset), channels, stride);
		}
		ConstView view(int offset=0) const {
			return ConstView(buffer.view(offset), channels, stride);
		}
		ConstView constView(int offset=0) const {
			return ConstView(buffer.view(offset), channels, stride);
		}

		MutableChannel operator[](int channel) {
			return buffer + channel*stride;
		}
		ConstChannel operator[](int channel) const {
			return buffer + channel*stride;
		}
		
		MultiBuffer & operator ++() {
			++buffer;
			return *this;
		}
		MultiBuffer & operator +=(int i) {
			buffer += i;
			return *this;
		}
		MutableView operator ++(int) {
			return MutableView(buffer++, channels, stride);
		}
		MutableView operator +(int i) {
			return MutableView(buffer + i, channels, stride);
		}
		ConstView operator +(int i) const {
			return ConstView(buffer + i, channels, stride);
		}
		MultiBuffer & operator --() {
			--buffer;
			return *this;
		}
		MultiBuffer & operator -=(int i) {
			buffer -= i;
			return *this;
		}
		MutableView operator --(int) {
			return MutableView(buffer--, channels, stride);
		}
		MutableView operator -(int i) {
			return MutableView(buffer - i, channels, stride);
		}
		ConstView operator -(int i) const {
			return ConstView(buffer - i, channels, stride);
		}
	};
	
	/** \defgroup Interpolators Interpolators
		\ingroup PitchDelayProcessor
		@{ */
	/// Nearest-neighbour interpolator
	/// \diagram{delay-random-access-nearest.svg,aliasing and maximum amplitude/delay errors for different input frequencies}
	template<typename Sample>
	struct InterpolatorNearest {
		static constexpr int inputLength = 1;
		static constexpr Sample latency = -0.5; // Because we're truncating, which rounds down too often
	
		template<class Data>
		static Sample fractional(const Data &data, Sample) {
			return data[0];
		}
	};
	/// Linear interpolator
	/// \diagram{delay-random-access-linear.svg,aliasing and maximum amplitude/delay errors for different input frequencies}
	template<typename Sample>
	struct InterpolatorLinear {
		static constexpr int inputLength = 2;
		static constexpr int latency = 0;
	
		template<class Data>
		static Sample fractional(const Data &data, Sample fractional) {
			Sample a = data[0], b = data[1];
			return a + fractional*(b - a);
		}
	};
	/// Spline cubic interpolator
	/// \diagram{delay-random-access-cubic.svg,aliasing and maximum amplitude/delay errors for different input frequencies}
	template<typename Sample>
	struct InterpolatorCubic {
		static constexpr int inputLength = 4;
		static constexpr int latency = 1;
	
		template<class Data>
		static Sample fractional(const Data &data, Sample fractional) {
			// Cubic interpolation
			Sample a = data[0], b = data[1], c = data[2], d = data[3];
			Sample cbDiff = c - b;
			Sample k1 = (c - a)*0.5;
			Sample k3 = k1 + (d - b)*0.5 - cbDiff*2;
			Sample k2 = cbDiff - k3 - k1;
			return b + fractional*(k1 + fractional*(k2 + fractional*k3)); // 16 ops total, not including the indexing
		}
	};

	// Efficient Algorithms and Structures for Fractional PitchDelayProcessor Filtering Based on Lagrange Interpolation
	// Franck 2009 https://www.aes.org/e-lib/browse.cfm?elib=14647
	namespace _franck_impl {
		template<typename Sample, int n, int low, int high>
		struct ProductRange {
			using Array = std::array<Sample, (n + 1)>;
			static constexpr int mid = (low + high)/2;
			using Left = ProductRange<Sample, n, low, mid>;
			using Right = ProductRange<Sample, n, mid + 1, high>;

			Left left;
			Right right;

			const Sample total;
			ProductRange(Sample x) : left(x), right(x), total(left.total*right.total) {}

			template<class Data>
			Sample calculateResult(Sample extraFactor, const Data &data, const Array &invFactors) {
				return left.calculateResult(extraFactor*right.total, data, invFactors)
					+ right.calculateResult(extraFactor*left.total, data, invFactors);
			}
		};
		template<typename Sample, int n, int index>
		struct ProductRange<Sample, n, index, index> {
			using Array = std::array<Sample, (n + 1)>;

			const Sample total;
			ProductRange(Sample x) : total(x - index) {}

			template<class Data>
			Sample calculateResult(Sample extraFactor, const Data &data, const Array &invFactors) {
				return extraFactor*data[index]*invFactors[index];
			}
		};
	}
	/** Fixed-order Lagrange interpolation.
	\diagram{interpolator-LagrangeN.svg,aliasing and amplitude/delay errors for different sizes}
	*/
	template<typename Sample, int n>
	struct InterpolatorLagrangeN {
		static constexpr int inputLength = n + 1;
		static constexpr int latency = (n - 1)/2;

		using Array = std::array<Sample, (n + 1)>;
		Array invDivisors;

		InterpolatorLagrangeN() {
			for (int j = 0; j <= n; ++j) {
				double divisor = 1;
				for (int k = 0; k < j; ++k) divisor *= (j - k);
				for (int k = j + 1; k <= n; ++k) divisor *= (j - k);
				invDivisors[j] = 1/divisor;
			}
		}

		template<class Data>
		Sample fractional(const Data &data, Sample fractional) const {
			constexpr int mid = n/2;
			using Left = _franck_impl::ProductRange<Sample, n, 0, mid>;
			using Right = _franck_impl::ProductRange<Sample, n, mid + 1, n>;

			Sample x = fractional + latency;

			Left left(x);
			Right right(x);

			return left.calculateResult(right.total, data, invDivisors) + right.calculateResult(left.total, data, invDivisors);
		}
	};
	template<typename Sample>
	using InterpolatorLagrange3 = InterpolatorLagrangeN<Sample, 3>;
	template<typename Sample>
	using InterpolatorLagrange7 = InterpolatorLagrangeN<Sample, 7>;
	template<typename Sample>
	using InterpolatorLagrange19 = InterpolatorLagrangeN<Sample, 19>;

	/** Fixed-size Kaiser-windowed sinc interpolation.
	\diagram{interpolator-KaiserSincN.svg,aliasing and amplitude/delay errors for different sizes}
	If `minimumPhase` is enabled, a minimum-phase version of the kernel is used:
	\diagram{interpolator-KaiserSincN-min.svg,aliasing and amplitude/delay errors for minimum-phase mode}
	*/
	template<typename Sample, int n, bool minimumPhase=false>
	struct InterpolatorKaiserSincN {
		static constexpr int inputLength = n;
		static constexpr Sample latency = minimumPhase ? 0 : (n*Sample(0.5) - 1);

		int subSampleSteps;
		std::vector<Sample> coefficients;
		
		InterpolatorKaiserSincN() : InterpolatorKaiserSincN(0.5 - 0.45/std::sqrt(n)) {}
		InterpolatorKaiserSincN(double passFreq) : InterpolatorKaiserSincN(passFreq, 1 - passFreq) {}
		InterpolatorKaiserSincN(double passFreq, double stopFreq) {
			subSampleSteps = 2*n; // Heuristic again.  Really it depends on the bandwidth as well.
			double kaiserBandwidth = (stopFreq - passFreq)*(n + 1.0/subSampleSteps);
			kaiserBandwidth += 1.25/kaiserBandwidth; // We want to place the first zero, but (because using this to window a sinc essentially integrates it in the freq-domain), our ripples (and therefore zeroes) are out of phase.  This is a heuristic fix.
			double sincScale = M_PI*(passFreq + stopFreq);

			double centreIndex = n*subSampleSteps*0.5, scaleFactor = 1.0/subSampleSteps;
			std::vector<Sample> windowedSinc(subSampleSteps*n + 1);
			
			::signalsmith::windows::Kaiser::withBandwidth(kaiserBandwidth, false).fill(windowedSinc, windowedSinc.size());

			for (size_t i = 0; i < windowedSinc.size(); ++i) {
				double x = (i - centreIndex)*scaleFactor;
				int intX = std::round(x);
				if (intX != 0 && std::abs(x - intX) < 1e-6) {
					// Exact 0s
					windowedSinc[i] = 0;
				} else if (std::abs(x) > 1e-6) {
					double p = x*sincScale;
					windowedSinc[i] *= std::sin(p)/p;
				}
			}
			
			if (minimumPhase) {
				signalsmith::fft::FFT<Sample> fft(windowedSinc.size()*2, 1);
				windowedSinc.resize(fft.size(), 0);
				std::vector<std::complex<Sample>> spectrum(fft.size());
				std::vector<std::complex<Sample>> cepstrum(fft.size());
				fft.fft(windowedSinc, spectrum);
				for (size_t i = 0; i < fft.size(); ++i) {
					spectrum[i] = std::log(std::abs(spectrum[i]) + 1e-30);
				}
				fft.fft(spectrum, cepstrum);
				for (size_t i = 1; i < fft.size()/2; ++i) {
					cepstrum[i] *= 0;
				}
				for (size_t i = fft.size()/2 + 1; i < fft.size(); ++i) {
					cepstrum[i] *= 2;
				}
				Sample scaling = Sample(1)/fft.size();
				fft.ifft(cepstrum, spectrum);

				for (size_t i = 0; i < fft.size(); ++i) {
					Sample phase = spectrum[i].imag()*scaling;
					Sample mag = std::exp(spectrum[i].real()*scaling);
					spectrum[i] = {mag*std::cos(phase), mag*std::sin(phase)};
				}
				fft.ifft(spectrum, cepstrum);
				windowedSinc.resize(subSampleSteps*n + 1);
				windowedSinc.shrink_to_fit();
				for (size_t i = 0; i < windowedSinc.size(); ++i) {
					windowedSinc[i] = cepstrum[i].real()*scaling;
				}
			}
			
			// Re-order into FIR fractional-delay blocks
			coefficients.resize(n*(subSampleSteps + 1));
			for (int k = 0; k <= subSampleSteps; ++k) {
				for (int i = 0; i < n; ++i) {
					coefficients[k*n + i] = windowedSinc[(subSampleSteps - k) + i*subSampleSteps];
				}
			}
		}
		
		template<class Data>
		Sample fractional(const Data &data, Sample fractional) const {
			Sample subSampleDelay = fractional*subSampleSteps;
			int lowIndex = subSampleDelay;
			if (lowIndex >= subSampleSteps) lowIndex = subSampleSteps - 1;
			Sample subSampleFractional = subSampleDelay - lowIndex;
			int highIndex = lowIndex + 1;
			
			Sample sumLow = 0, sumHigh = 0;
			const Sample *coeffLow = coefficients.data() + lowIndex*n;
			const Sample *coeffHigh = coefficients.data() + highIndex*n;
			for (int i = 0; i < n; ++i) {
				sumLow += data[i]*coeffLow[i];
				sumHigh += data[i]*coeffHigh[i];
			}
			return sumLow + (sumHigh - sumLow)*subSampleFractional;
		}
	};

	template<typename Sample>
	using InterpolatorKaiserSinc20 = InterpolatorKaiserSincN<Sample, 20>;
	template<typename Sample>
	using InterpolatorKaiserSinc8 = InterpolatorKaiserSincN<Sample, 8>;
	template<typename Sample>
	using InterpolatorKaiserSinc4 = InterpolatorKaiserSincN<Sample, 4>;

	template<typename Sample>
	using InterpolatorKaiserSinc20Min = InterpolatorKaiserSincN<Sample, 20, true>;
	template<typename Sample>
	using InterpolatorKaiserSinc8Min = InterpolatorKaiserSincN<Sample, 8, true>;
	template<typename Sample>
	using InterpolatorKaiserSinc4Min = InterpolatorKaiserSincN<Sample, 4, true>;
	///  @}
	
	/** @brief A delay-line reader which uses an external buffer
 
		This is useful if you have multiple delay-lines reading from the same buffer.
	*/
	template<class Sample, template<typename> class Interpolator=InterpolatorLinear>
	class Reader : public Interpolator<Sample> /* so we can get the empty-base-class optimisation */ {
		using Super = Interpolator<Sample>;
	public:
		Reader () {}
		/// Pass in a configured interpolator
		Reader (const Interpolator<Sample> &interpolator) : Super(interpolator) {}
	
		template<typename Buffer>
		Sample read(const Buffer &buffer, Sample delaySamples) const {
			int startIndex = delaySamples;
			Sample remainder = delaySamples - startIndex;
			
			// PitchDelayProcessor buffers use negative indices, but interpolators use positive ones
			using View = decltype(buffer - startIndex);
			struct Flipped {
				 View view;
				 Sample operator [](int i) const {
					return view[-i];
				 }
			};
			return Super::fractional(Flipped{buffer - startIndex}, remainder);
		}
	};

	/**	@brief A single-channel delay-line containing its own buffer.*/
	template<class Sample, template<typename> class Interpolator=InterpolatorLinear>
	class Delay : private Reader<Sample, Interpolator> {
		using Super = Reader<Sample, Interpolator>;
		Buffer<Sample> buffer;
	public:
		static constexpr Sample latency = Super::latency;

		Delay(int capacity=0) : buffer(1 + capacity + Super::inputLength) {}
		/// Pass in a configured interpolator
		Delay(const Interpolator<Sample> &interp, int capacity=0) : Super(interp), buffer(1 + capacity + Super::inputLength) {}
		
		void reset(Sample value=Sample()) {
			buffer.reset(value);
		}
		void resize(int minCapacity, Sample value=Sample()) {
			buffer.resize(minCapacity + Super::inputLength, value);
		}
		
		/** Read a sample from `delaySamples` >= 0 in the past.
		The interpolator may add its own latency on top of this (see `PitchDelayProcessor::latency`).  The default interpolation (linear) has 0 latency.
		*/
		Sample read(Sample delaySamples) const {
			return Super::read(buffer, delaySamples);
		}
		/// Writes a sample. Returns the same object, so that you can say `delay.write(v).read(delay)`.
		Delay & write(Sample value) {
			++buffer;
			buffer[0] = value;
			return *this;
		}
	};

	/**	@brief A multi-channel delay-line with its own buffer. */
	template<class Sample, template<typename> class Interpolator=InterpolatorLinear>
	class MultiDelay : private Reader<Sample, Interpolator> {
		using Super = Reader<Sample, Interpolator>;
		int channels;
		MultiBuffer<Sample> multiBuffer;
	public:
		static constexpr Sample latency = Super::latency;

		MultiDelay(int channels=0, int capacity=0) : channels(channels), multiBuffer(channels, 1 + capacity + Super::inputLength) {}

		void reset(Sample value=Sample()) {
			multiBuffer.reset(value);
		}
		void resize(int nChannels, int capacity, Sample value=Sample()) {
			channels = nChannels;
			multiBuffer.resize(channels, capacity + Super::inputLength, value);
		}
		
		/// A single-channel delay-line view, similar to a `const PitchDelayProcessor`
		struct ChannelView {
			static constexpr Sample latency = Super::latency;

			const Super &reader;
			typename MultiBuffer<Sample>::ConstChannel channel;
			
			Sample read(Sample delaySamples) const {
				return reader.read(channel, delaySamples);
			}
		};
		ChannelView operator [](int channel) const {
			return ChannelView{*this, multiBuffer[channel]};
		}

		/// A multi-channel result, lazily calculating samples
		struct DelayView {
			Super &reader;
			typename MultiBuffer<Sample>::ConstView view;
			Sample delaySamples;
			
			// Calculate samples on-the-fly
			Sample operator [](int c) const {
				return reader.read(view[c], delaySamples);
			}
		};
		DelayView read(Sample delaySamples) {
			return DelayView{*this, multiBuffer.constView(), delaySamples};
		}
		/// Reads into the provided output structure
		template<class Output>
		void read(Sample delaySamples, Output &output) {
			for (int c = 0; c < channels; ++c) {
				output[c] = Super::read(multiBuffer[c], delaySamples);
			}
		}
		/// Reads separate delays for each channel
		template<class Delays, class Output>
		void readMulti(const Delays &delays, Output &output) {
			for (int c = 0; c < channels; ++c) {
				output[c] = Super::read(multiBuffer[c], delays[c]);
			}
		}
		template<class Data>
		MultiDelay & write(const Data &data) {
			++multiBuffer;
			for (int c = 0; c < channels; ++c) {
				multiBuffer[c][0] = data[c];
			}
			return *this;
		}
	};

/** @} */
}
    // spectral.h
    namespace spectral {
	/**	@defgroup Spectral Spectral Processing
		@brief Tools for frequency-domain manipulation of audio signals
		
		@{
		@file
	*/
	
	/** @brief An FFT with built-in windowing and round-trip scaling
	
		This uses a Modified Real FFT, which applies half-bin shift before the transform.  The result therefore has `N/2` bins, centred at the frequencies: `(i + 0.5)/N`.
		
		This avoids the awkward (real-valued) bands for DC-offset and Nyquist.
	 */
	template<typename Sample>
	class WindowedFFT {
		using MRFFT = signalsmith::fft::ModifiedRealFFT<Sample>;
		using Complex = std::complex<Sample>;
		MRFFT mrfft{2};

		std::vector<Sample> fftWindow;
		std::vector<Sample> timeBuffer;
		std::vector<Complex> freqBuffer;
	public:
		/// Returns a fast FFT size <= `size`
		static int fastSizeAbove(int size, int divisor=1) {
			return MRFFT::fastSizeAbove(size/divisor)*divisor;
		}
		/// Returns a fast FFT size >= `size`
		static int fastSizeBelow(int size, int divisor=1) {
			return MRFFT::fastSizeBelow(1 + (size - 1)/divisor)*divisor;
		}

		WindowedFFT() {}
		WindowedFFT(int size) {
			setSize(size);
		}
		template<class WindowFn>
		WindowedFFT(int size, WindowFn fn, Sample windowOffset=0.5) {
			setSize(size, fn, windowOffset);
		}

		/// Sets the size, returning the window for modification (initially all 1s)
		std::vector<Sample> & setSizeWindow(int size) {
			mrfft.setSize(size);
			fftWindow.resize(size, 1);
			timeBuffer.resize(size);
			freqBuffer.resize(size);
			return fftWindow;
		}
		/// Sets the FFT size, with a user-defined functor for the window
		template<class WindowFn>
		void setSize(int size, WindowFn fn, Sample windowOffset=0.5) {
			setSizeWindow(size);
		
			Sample invSize = 1/(Sample)size;
			for (int i = 0; i < size; ++i) {
				Sample r = (i + windowOffset)*invSize;
				fftWindow[i] = fn(r);
			}
		}
		/// Sets the size (using the default Blackman-Harris window)
		void setSize(int size) {
			setSize(size, [](double x) {
				double phase = 2*M_PI*x;
				// Blackman-Harris
				return 0.35875 + 0.48829*std::cos(phase) + 0.14128*std::cos(phase*2) + 0.1168*std::cos(phase*3);
			});
		}

		const std::vector<Sample> & window() const {
			return this->fftWindow;
		}
		int size() const {
			return mrfft.size();
		}
		
		/// Performs an FFT (with windowing)
		template<class Input, class Output>
		void fft(Input &&input, Output &&output) {
			struct WindowedInput {
				const Input &input;
				std::vector<Sample> &window;
				SIGNALSMITH_INLINE Sample operator [](int i) {
					return input[i]*window[i];
				}
			};
		
			mrfft.fft(WindowedInput{input, fftWindow}, output);
		}
		/// Performs an FFT (no windowing)
		template<class Input, class Output>
		void fftRaw(Input &&input, Output &&output) {
			mrfft.fft(input, output);
		}

		/// Inverse FFT, with windowing and 1/N scaling
		template<class Input, class Output>
		void ifft(Input &&input, Output &&output) {
			mrfft.ifft(input, output);
			int size = mrfft.size();
			Sample norm = 1/(Sample)size;
			for (int i = 0; i < size; ++i) {
				output[i] *= norm*fftWindow[i];
			}
		}
		/// Performs an IFFT (no windowing)
		template<class Input, class Output>
		void ifftRaw(Input &&input, Output &&output) {
			mrfft.ifft(input, output);
		}
	};
	
	/** STFT synthesis, built on a `MultiBuffer`.
 
		Any window length and block interval is supported, but the FFT size may be rounded up to a faster size (by zero-padding).  It uses a heuristically-optimal Kaiser window modified for perfect-reconstruction.
		
		\diagram{stft-aliasing-simulated.svg,Simulated bad-case aliasing (random phase-shift for each band) for overlapping ratios}

		There is a "latest valid index", and you can read the output up to one `historyLength` behind this (see `.resize()`).  You can read up to one window-length _ahead_ to get partially-summed future output.
		
		\diagram{stft-buffer-validity.svg}
		
		You move the valid index along using `.ensureValid()`, passing in a functor which provides spectra (using `.analyse()` and/or direct modification through `.spectrum[c]`):

		\code
			void processSample(...) {
				stft.ensureValid([&](int) {
					// Here, we introduce (1 - windowSize) of latency
					stft.analyse(inputBuffer.view(1 - windowSize))
				});
				// read as a MultiBuffer
				auto result = stft.at(0);
				++stft; // also moves the latest valid index
			}

			void processBlock(...) {
				// assuming `historyLength` == blockSize
				stft.ensureValid(blockSize, [&](int blockStartIndex) {
					int inputStart = blockStartIndex + (1 - windowSize);
					stft.analyse(inputBuffer.view(inputStart));
				});
				auto earliestValid = stft.at(0);
				auto latestValid = stft.at(blockSize);
				stft += blockSize;
			}
		\endcode
		
		The index passed to this functor will be greater than the previous valid index, and `<=` the index you pass in.  Therefore, if you call `.ensureValid()` every sample, it can only ever be `0`.
	*/
	template<typename Sample>
	class STFT : public signalsmith::delay::MultiBuffer<Sample> {
		using Super = signalsmith::delay::MultiBuffer<Sample>;
		using Complex = std::complex<Sample>;

		int channels = 0, _windowSize = 0, _fftSize = 0, _interval = 1;
		int validUntilIndex = 0;

		class MultiSpectrum {
			int channels, stride;
			std::vector<Complex> buffer;
		public:
			MultiSpectrum() : MultiSpectrum(0, 0) {}
			MultiSpectrum(int channels, int bands) : channels(channels), stride(bands), buffer(channels*bands, 0) {}
			
			void resize(int nChannels, int nBands) {
				channels = nChannels;
				stride = nBands;
				buffer.assign(channels*stride, 0);
			}
			
			void reset() {
				buffer.assign(buffer.size(), 0);
			}
			
			void swap(MultiSpectrum &other) {
				using std::swap;
				swap(buffer, other.buffer);
			}

			Complex * operator [](int channel) {
				return buffer.data() + channel*stride;
			}
			const Complex * operator [](int channel) const {
				return buffer.data() + channel*stride;
			}

            [[maybe_unused]] [[nodiscard]] SIGNALSMITH_INLINE int getChannels() const { return channels; }
            [[maybe_unused]] [[nodiscard]] SIGNALSMITH_INLINE int getStride() const { return stride; }
		};
		std::vector<Sample> timeBuffer;

		void resizeInternal(int newChannels, int windowSize, int newInterval, int historyLength, int zeroPadding) {
			Super::resize(newChannels,
				windowSize /* for output summing */
				+ newInterval /* so we can read `windowSize` ahead (we'll be at most `interval-1` from the most recent block */
				+ historyLength);

			int fftSize = fft.fastSizeAbove(windowSize + zeroPadding);
			
			this->channels = newChannels;
			_windowSize = windowSize;
			this->_fftSize = fftSize;
			this->_interval = newInterval;
			validUntilIndex = -1;

			auto &window = fft.setSizeWindow(fftSize);
			if (windowShape == Window::kaiser) {
				using Kaiser = ::signalsmith::windows::Kaiser;
				/// Roughly optimal Kaiser for STFT analysis (forced to perfect reconstruction)
				auto kaiser = Kaiser::withBandwidth(windowSize/double(_interval), true);
				kaiser.fill(window, windowSize);
			} else {
				using Confined = ::signalsmith::windows::ApproximateConfinedGaussian;
				auto confined = Confined::withBandwidth(windowSize/double(_interval));
				confined.fill(window, windowSize);
			}
			::signalsmith::windows::forcePerfectReconstruction(window, windowSize, _interval);
			
			// TODO: fill extra bits of an input buffer with NaN/Infinity, to break this, and then fix by adding zero-padding to WindowedFFT (as opposed to zero-valued window sections)
			for (int i = windowSize; i < fftSize; ++i) {
				window[i] = 0;
			}

			spectrum.resize(channels, fftSize/2);
			timeBuffer.resize(fftSize);
		}
	public:
		/** Swaps between the default (Kaiser) shape and Approximate Confined Gaussian (ACG).
		\diagram{stft-windows.svg,Default (Kaiser) windows and partial cumulative sum}
		The ACG has better rolloff since its edges go to 0:
		\diagram{stft-windows-acg.svg,ACG windows and partial cumulative sum}
		However, it generally has worse performance in terms of total sidelobe energy, affecting worst-case aliasing levels for (most) higher overlap ratios:
		\diagram{stft-aliasing-simulated-acg.svg,Simulated bad-case aliasing for ACG windows - compare with above}*/
		enum class Window {kaiser, acg};
		Window windowShape = Window::kaiser;
		
		using Spectrum = MultiSpectrum;
		Spectrum spectrum;
		WindowedFFT<Sample> fft;
		
		STFT() {}
		/// Parameters passed straight to `.resize()`
		STFT(int channels, int windowSize, int interval, int historyLength=0, int zeroPadding=0) {
			resize(channels, windowSize, interval, historyLength, zeroPadding);
		}

		/// Sets the channel-count, FFT size and interval.
		void resize(int nChannels, int windowSize, int interval, int historyLength=0, int zeroPadding=0) {
			resizeInternal(nChannels, windowSize, interval, historyLength, zeroPadding);
		}
		
		int windowSize() const {
			return _windowSize;
		}
		int fftSize() const {
			return _fftSize;
		}
		int interval() const {
			return _interval;
		}
		/// Returns the (analysis and synthesis) window
		decltype(fft.window()) window() const {
			return fft.window();
		}
		/// Calculates the effective window for the partially-summed future output (relative to the most recent block)
		std::vector<Sample> partialSumWindow() const {
			const auto &w = window();
			std::vector<Sample> result(_windowSize, 0);
			for (int offset = 0; offset < _windowSize; offset += _interval) {
				for (int i = 0; i < _windowSize - offset; ++i) {
					Sample value = w[i + offset];
					result[i] += value*value;
				}
			}
			return result;
		}
		
		/// Resets everything - since we clear the output sum, it will take `windowSize` samples to get proper output.
		void reset() {
			Super::reset();
			spectrum.reset();
			validUntilIndex = -1;
		}
		
		/** Generates valid output up to the specified index (or 0), using the callback as many times as needed.
			
			The callback should be a functor accepting a single integer argument, which is the index for which a spectrum is required.
			
			The block created from these spectra will start at this index in the output, plus `.latency()`.
		*/
		template<class AnalysisFn>
		void ensureValid(int i, AnalysisFn fn) {
			while (validUntilIndex < i) {
				int blockIndex = validUntilIndex + 1;
				fn(blockIndex);

				auto output = this->view(blockIndex);
				for (int c = 0; c < channels; ++c) {
					auto channel = output[c];

					// Clear out the future sum, a window-length and an interval ahead
					for (int wi = _windowSize; wi < _windowSize + _interval; ++wi) {
						channel[wi] = 0;
					}

					// Add in the IFFT'd result
					fft.ifft(spectrum[c], timeBuffer);
					for (int wi = 0; wi < _windowSize; ++wi) {
						channel[wi] += timeBuffer[wi];
					}
				}
				validUntilIndex += _interval;
			}
		}
		/// The same as above, assuming index 0
		template<class AnalysisFn>
		void ensureValid(AnalysisFn fn) {
			return ensureValid(0, fn);
		}
		
		/** Analyse a multi-channel input, for any type where `data[channel][index]` returns samples
 
		Results can be read/edited using `.spectrum`. */
		template<class Data>
		void analyse(Data &&data) {
			for (int c = 0; c < channels; ++c) {
				fft.fft(data[c], spectrum[c]);
			}
		}
		template<class Data>
		void analyse(int c, Data &&data) {
			fft.fft(data, spectrum[c]);
		}
		/// Analyse without windowing
		template<class Data>
		void analyseRaw(Data &&data) {
			for (int c = 0; c < channels; ++c) {
				fft.fftRaw(data[c], spectrum[c]);
			}
		}
		template<class Data>
		void analyseRaw(int c, Data &&data) {
			fft.fftRaw(data, spectrum[c]);
		}

		int bands() const {
			return _fftSize/2;
		}

		/** Internal latency (between the block-index requested in `.ensureValid()` and its position in the output)
 
		Currently unused, but it's in here to allow for a future implementation which spreads the FFT calculations out across each interval.*/
		int latency() {
			return 0;
		}
		
		// @name Shift the underlying buffer (moving the "valid" index accordingly)
		// @{
		STFT & operator ++() {
			Super::operator ++();
			validUntilIndex--;
			return *this;
		}
		STFT & operator +=(int i) {
			Super::operator +=(i);
			validUntilIndex -= i;
			return *this;
		}
		STFT & operator --() {
			Super::operator --();
			validUntilIndex++;
			return *this;
		}
		STFT & operator -=(int i) {
			Super::operator -=(i);
			validUntilIndex += i;
			return *this;
		}
		// @}

		typename Super::MutableView operator ++(int postIncrement) {
			auto result = Super::operator ++(postIncrement);
			validUntilIndex--;
			return result;
		}
		typename Super::MutableView operator --(int postIncrement) {
			auto result = Super::operator --(postIncrement);
			validUntilIndex++;
			return result;
		}
	};

	/** STFT processing, with input/output.
		Before calling `.ensureValid(index)`, you should make sure the input is filled up to `index`.
	*/
	template<typename Sample>
	class ProcessSTFT : public STFT<Sample> {
		using Super = STFT<Sample>;
	public:
		signalsmith::delay::MultiBuffer<Sample> input;
	
		ProcessSTFT(int inChannels, int outChannels, int windowSize, int interval, int historyLength=0) {
			resize(inChannels, outChannels, windowSize, interval, historyLength);
		}

		/** Alter the spectrum, using input up to this point, for the output block starting from this point.
			Sub-classes should replace this with whatever processing is desired. */
		virtual void processSpectrum(int /*blockIndex*/) {}
		
		/// Sets the input/output channels, FFT size and interval.
		void resize(int inChannels, int outChannels, int windowSize, int interval, int historyLength=0) {
			Super::resize(outChannels, windowSize, interval, historyLength);
			input.resize(inChannels, windowSize + interval + historyLength);
		}
		void reset(Sample value=Sample()) {
			Super::reset(value);
			input.reset(value);
		}

		/// Internal latency, including buffering samples for analysis.
		int latency() {
			return Super::latency() + (this->windowSize() - 1);
		}
		
		void ensureValid(int i=0) {
			Super::ensureValid(i, [&](int blockIndex) {
				this->analyse(input.view(blockIndex - this->windowSize() + 1));
				this->processSpectrum(blockIndex);
			});
		}

		// @name Shift the output, input, and valid index.
		// @{
		ProcessSTFT & operator ++() {
			Super::operator ++();
			++input;
			return *this;
		}
		ProcessSTFT & operator +=(int i) {
			Super::operator +=(i);
			input += i;
			return *this;
		}
		ProcessSTFT & operator --() {
			Super::operator --();
			--input;
			return *this;
		}
		ProcessSTFT & operator -=(int i) {
			Super::operator -=(i);
			input -= i;
			return *this;
		}
		// @}
	};

/** @} */
}
} 