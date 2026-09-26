#pragma once

#include <charconv>
#include <string_view>
#include <system_error>

//Parses all of str as a base-10 number, ignoring any spaces or tabs around it.
//Returns false and leaves value unchanged if str holds anything else or the number does not fit
//in T.
template<typename T>
bool parse_number(std::string_view str, T& value) {
	while (!str.empty() && (str.front() == ' ' || str.front() == '\t')) str.remove_prefix(1);
	while (!str.empty() && (str.back() == ' ' || str.back() == '\t')) str.remove_suffix(1);
	const char* first = str.data();
	const char* last = first + str.size();
	T result{};
	auto [ptr, ec] = std::from_chars(first, last, result);
	if (ec != std::errc() || ptr != last) {
		return false;
	}
	value = result;
	return true;
}
