#include <QApplication>
#include <iostream>
#include <QLocale>
#include <QTranslator>
#include <QDebug>
#include <QElapsedTimer>
#include <QSplashScreen>
#include <QSurfaceFormat>

#include "../libs/basic/logger.h"
#include "main_window.h"

#include <time.h>



class Application : public QApplication
{
public:
	Application(int& argc, char ** argv) : QApplication(argc, argv) { }
	virtual ~Application() { }

	// reimplemented from QApplication so we can throw exceptions in slots
	virtual bool notify(QObject * receiver, QEvent * event) {
		try {
			return QApplication::notify(receiver, event);
		}
		catch (std::exception& e) {
			qCritical() << "Exception thrown:" << e.what();
		}
		return false;
	}
};


int main(int argc, char **argv)
{
    srand(time(nullptr));

    //Locale management
    //Force 'English' locale to get a consistent behavior everywhere
    QLocale locale = QLocale(QLocale::English);
    locale.setNumberOptions(QLocale::c().numberOptions());
    QLocale::setDefault(locale);

#ifdef Q_OS_UNIX
    //We reset the numeric locale for POSIX functions
    //See http://qt-project.org/doc/qt-5/qcoreapplication.html#locale-settings
    setlocale(LC_NUMERIC, "C");
#endif

#if (QT_VERSION >= QT_VERSION_CHECK(5, 6, 0) && (QT_VERSION < QT_VERSION_CHECK(6, 0, 0)))
    QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
#endif

    // Note: Calling QSurfaceFormat::setDefaultFormat() before constructing the
    //       QApplication instance is mandatory on some platforms(for example, macOS)
    //       when an OpenGL core profile context is requested. This is to ensure
    //       that resource sharing between contexts stays functional as all internal
    //       contexts are created using the correct version and profile.
    QSurfaceFormat format = QSurfaceFormat::defaultFormat();
    format.setVersion(4, 1);
    format.setProfile(QSurfaceFormat::CompatibilityProfile);
    format.setDepthBufferSize(24);
    format.setStencilBufferSize(8);
    format.setSwapBehavior(QSurfaceFormat::DoubleBuffer);
    format.setSamples(4);
#ifndef NDEBUG
    format.setOption(QSurfaceFormat::DebugContext);
#endif
    QSurfaceFormat::setDefaultFormat(format);

    Application app(argc, argv);

#ifndef _DEBUG
	// splash screen
	QPixmap pixmap(":/Resources/splash.png");
	QSplashScreen splash(pixmap, Qt::WindowStaysOnTopHint);
    QElapsedTimer splashTimer;
	splashTimer.start();
	splash.show();
	splash.showMessage("  Starting MVStudio...");
    QApplication::processEvents();

	//we want the splash screen to be visible a minimum amount of time (in milliseconds.)
	while (splashTimer.elapsed() < 500) {
		splash.raise();
		QApplication::processEvents(); //to let the system breath!
	}
#endif

	MainWindow window;	
	window.show();

#ifndef _DEBUG
	splash.finish(&window);
	//splash.close();
	QApplication::processEvents();
#endif

	//if (VerifyRegistry() == false)
	//	return 1;

	return app.exec();
};
