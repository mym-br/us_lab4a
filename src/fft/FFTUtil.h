#ifndef FFTUTIL_H
#define FFTUTIL_H

#include "Exception.h"



namespace Lab {
namespace FFTUtil {

unsigned int nextFastSize(unsigned int n);
unsigned int nextFastEvenSize(unsigned int n);



inline
unsigned int
nextFastSize(unsigned int n)
{
	// Copied from Kiss FFT and slightly modified.
	/*
	Copyright (c) 2003-2010, Mark Borgerding

	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions
	are met:

	* Redistributions of source code must retain the above copyright
	  notice, this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright
	  notice, this list of conditions and the following disclaimer in the
	  documentation and/or other materials provided with the distribution.
	* Neither the author nor the names of any contributors may be used to
	  endorse or promote products derived from this software without
	  specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
	OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
	LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
	THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	*/
	if (n == 0) {
		THROW_EXCEPTION(InvalidParameterException, "n must be > 0.");
	}
	if (n == 1) n = 2;
	while (1) {
		unsigned int m = n;
		while ((m & 1) == 0) m >>= 1;
		while (m % 3 == 0) m /= 3;
		while (m % 5 == 0) m /= 5;
		if (m == 1) {
			break; // n is completely factorable by twos, threes, and fives
		}
		++n;
	}
	return n;
}

inline
unsigned int
nextFastEvenSize(unsigned int n)
{
	// Copied from Kiss FFT and modified.
	/*
	Copyright (c) 2003-2010, Mark Borgerding

	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions
	are met:

	* Redistributions of source code must retain the above copyright
	  notice, this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright
	  notice, this list of conditions and the following disclaimer in the
	  documentation and/or other materials provided with the distribution.
	* Neither the author nor the names of any contributors may be used to
	  endorse or promote products derived from this software without
	  specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
	OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
	LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
	THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	*/
	if (n == 0) {
		THROW_EXCEPTION(InvalidParameterException, "n must be > 0.");
	}
	if (n & 1) ++n;
	while (1) {
		unsigned int m = n;
		while ((m & 1) == 0) m >>= 1;
		while (m % 3 == 0) m /= 3;
		while (m % 5 == 0) m /= 5;
		if (m == 1) {
			break; // n is completely factorable by twos, threes, and fives and is even
		}
		n += 2;
	}
	return n;
}

} // namespace FFTUtil
} // namescape Lab







#endif // FFTUTIL_H
