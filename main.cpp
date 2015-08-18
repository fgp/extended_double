#include <stdlib.h>
#include <stdexcept>

#include <boost/test/unit_test.hpp>
#include <boost/test/prg_exec_monitor.hpp>

int test_main(int, char*[]);
int cpp_main(int argc, char** argv);

/**
 * \fn init_unit_test_suite
 *
 * Unit test initialization function. Supplied by boost::unit_test.
 *
 * Declaration copied from `boost/test/impl/unit_test_main.ipp`
 */
extern ::boost::unit_test::test_suite* init_unit_test_suite( int argc, char* argv[] );

/**
 * \fn test_main
 *
 * For wheatver reason, boost::unit_test thinks the use will supply this, and
 * linking fails if it isn't present. A suitable combination of BOOST_TEST_MAIN,
 * BOOST_TEST_NOMAIN and BOOST_TEST_MODULE might be able to convince boost
 * not to assume this function is present, but I gave up trying to understand
 * how exactly all these flags interact, and decided that simply supplying the
 * function is easier.
 *
 * Declaration copied from `boost/test/impl/test_main.ipp`
 */
int test_main(int, char*[]) {
	return 0;
}

/**
 * \fn cpp_main
 *
 * Copied from `boost/test/impl/unit_test_main.ipp`
 */
int cpp_main(int argc, char** argv) {
	boost::unit_test::init_unit_test_func init_func = &init_unit_test_suite;
	return ::boost::unit_test::unit_test_main(init_func, argc, argv);
}

/**
 * \fn main
 *
 * main() function that calls the boost::unit_test test runner. It should
 * be possible to make boost supply this method, but I don't know how, and
 * even if I did, I wouldn't know how to make it open CoverStory at the end
 * to view coverage results.
 *
 * Copied from `boost/test/prg_exec_monitor.hpp`
 */
int main(int argc, char** argv)
{
	const int r = ::boost::prg_exec_monitor_main(&cpp_main, argc, argv);

	return r;
}
