class Pair{
public:
	Pair(long i, long j):first(i), second(j){

	}
	long first;
	long second;
private:
	Pair(const Pair &other);
	Pair& operator=(const Pair& other);
};