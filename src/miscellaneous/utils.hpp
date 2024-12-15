#ifndef UTILS_HPP
#define UTILS_HPP 1
/*!
\brief Returns a vector containing only unique elements
\param v The input vector, must be SORTED!!
\return The unique vector
*/
template <class T> vector<T> unique(vector<T> const& v);

template <class T> vector<T> unique(vector<T> const& v)
{
	std::size_t n = v.size();
	vector<T> res;
	res.reserve(n);
	if (n == 0)
		return res;
	res.push_back(v[0]);
	for (typename vector<T>::const_iterator it = v.begin() + 1; it != v.end(); ++it)
		if (*it == *(it - 1))
			continue;
		else
			res.push_back(*it);
	return res;
}

#endif // UTILS_HPP