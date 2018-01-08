#ifndef MINSTDPSEUDORANDOMNUMBERGENERATOR_H_
#define MINSTDPSEUDORANDOMNUMBERGENERATOR_H_



namespace Lab {

class MinstdPseudorandomNumberGenerator {
public:
	MinstdPseudorandomNumberGenerator(long seed);
	~MinstdPseudorandomNumberGenerator();

	// Returns one value of the pseudorandom sequence ]0.0,1.0[.
	double get();
private:
	//MinstdPseudorandomNumberGenerator(const MinstdPseudorandomNumberGenerator&);
	//MinstdPseudorandomNumberGenerator& operator=(const MinstdPseudorandomNumberGenerator&);

	static const double a;
	static const double m;
	static const double invM;
	double x_;
};

} // namespace Lab

#endif /* MINSTDPSEUDORANDOMNUMBERGENERATOR_H_ */
